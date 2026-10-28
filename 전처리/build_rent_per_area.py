#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_rent_per_area.py
======================

영국 우편번호 앞자리(outcode, 예: BA2) 단위로 "단위면적당 평균 임대료"를 산출한다.

출력 CSV 컬럼 (필수 두 개 포함):
    outcode             우편번호 앞자리 (예: BA2)        <-- 필수
    rent_per_sqm        단위면적당 평균 월 임대료 (£/㎡)  <-- 필수
    rent_per_sqft       참고용 (£/sqft)
    rent_per_pyeong     참고용 (£/평, 1평=3.305785㎡)
    median_rent         해당 outcode에 매핑된 LA의 중위 월세 (£)
    median_floor_area   해당 outcode의 중위 전용면적 (㎡, EPC 기준)
    local_authority     매핑된 지방자치단체 코드
    n_epc               면적 계산에 쓰인 EPC 표본 수
    low_sample          표본 부족 플래그 (n_epc < MIN_SAMPLE)
출력에는 영국의 (사용 가능한) 모든 outcode가 한 행씩 들어간다.
데이터가 없는 outcode도 행은 남고 임대료 값은 NaN 으로 표시된다.

----------------------------------------------------------------------
데이터 소스 (모두 무료, 직접 다운로드 필요)
----------------------------------------------------------------------
1) 면적 (EPC, 잉글랜드·웨일스)
   https://epc.opendatacommunities.org/  (무료 가입 후 'all results' bulk download)
   압축을 풀면 지방자치단체별 폴더 안에 certificates.csv 가 들어 있다.
   사용 컬럼: POSTCODE, TOTAL_FLOOR_AREA(㎡), LOCAL_AUTHORITY,
             NUMBER_HABITABLE_ROOMS, PROPERTY_TYPE
   * 스코틀랜드는 별도 등록부(scottishepcregister.org.uk) — 같은 형식으로 합치면 됨.

2) 임대료 (지방자치단체 단위)
   - 현재: ONS PIPR price levels
     https://www.ons.gov.uk/economy/inflationandpriceindices/bulletins/privaterentandhousepricesuk/latest
   - 과거 상세(LA × 침실수 중위 임대료): VOA PRMS (2023-09 까지, 아카이브)
     https://www.ons.gov.uk/.../privaterentalmarketsummarystatisticsinengland
   둘 중 무엇을 쓰든, 아래 형식의 '정규화된 CSV' 하나로 만들어서 --rent 로 넣는다:
       la_code, median_rent           (전체 기준)            ... 또는
       la_code, bedrooms, median_rent (침실수별 — 권장)
   * 정규화 방법은 README 참고 (PRMS 엑셀 → tidy CSV).

3) outcode 마스터 목록 (모든 outcode 보장용, 선택)
   ONS Postcode Directory (ONSPD), Open Geography Portal 무료 다운로드.
   사용 컬럼: pcds(전체 우편번호) — 여기서 outcode 추출.
   --onspd 를 주지 않으면 EPC 에 등장한 outcode 만 마스터로 사용한다.

----------------------------------------------------------------------
임대료 지리 단위에 대한 솔직한 한계
----------------------------------------------------------------------
공식 임대료 데이터의 최소 단위는 지방자치단체(LA)다. outcode 는 LA 보다 작고
하나의 LA 안에 여러 outcode 가 있으므로, 같은 LA 에 속한 outcode 들은 같은
median_rent 값을 공유한다. 따라서 최종 rent_per_sqm 의 outcode 간 차이는
주로 'outcode 별 전용면적 차이'에서 나온다. outcode 단위로 임대료 자체가
달라지길 원한다면 매물 스크래핑(rent+면적)을 augment 로 붙여야 한다
(--listings 옵션, 아래 참고).
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------
# 설정값
# ----------------------------------------------------------------------
SQM_PER_SQFT = 0.09290304          # 1 sqft = 0.092903 ㎡
SQM_PER_PYEONG = 3.3057851         # 1 평  = 3.305785 ㎡
MIN_SAMPLE = 30                    # outcode 면적 중위값 신뢰 최소 표본
FLOOR_AREA_MIN, FLOOR_AREA_MAX = 10.0, 1000.0   # ㎡, 이상치 컷오프

# EPC 원본 컬럼명 (잉글랜드·웨일스 표준 스키마)
EPC_COLS = {
    "postcode": "POSTCODE",
    "floor_area": "TOTAL_FLOOR_AREA",
    "la": "LOCAL_AUTHORITY",
    "rooms": "NUMBER_HABITABLE_ROOMS",
    "ptype": "PROPERTY_TYPE",
}

OUTCODE_RE = re.compile(r"^([A-Z]{1,2}\d[A-Z\d]?)\s*\d[A-Z]{2}$")


# ----------------------------------------------------------------------
# 유틸: 우편번호 -> outcode (앞자리)
# ----------------------------------------------------------------------
def to_outcode(postcode: object) -> str | float:
    """'BA2 7AY' / 'ba27ay' -> 'BA2'. 유효하지 않으면 NaN."""
    if not isinstance(postcode, str):
        return np.nan
    pc = postcode.upper().strip().replace("  ", " ")
    pc_nospace = pc.replace(" ", "")
    if len(pc_nospace) < 5:
        return np.nan
    # 표준 우편번호는 항상 마지막 3자가 inward code(숫자+영문2)
    candidate = pc_nospace[:-3]
    m = re.match(r"^[A-Z]{1,2}\d[A-Z\d]?$", candidate)
    return candidate if m else np.nan


# ----------------------------------------------------------------------
# 1) EPC 로딩 -> outcode 별 중위 면적 + outcode->LA 매핑
# ----------------------------------------------------------------------
def iter_epc_files(epc_path: str):
    """디렉터리면 **/certificates.csv 들을, 단일 CSV면 그 파일을 순회."""
    if os.path.isdir(epc_path):
        files = sorted(glob.glob(os.path.join(epc_path, "**", "certificates.csv"),
                                 recursive=True))
        if not files:                       # 평범한 csv 모음일 수도
            files = sorted(glob.glob(os.path.join(epc_path, "**", "*.csv"),
                                     recursive=True))
        for f in files:
            yield f
    else:
        yield epc_path


def load_epc(epc_path: str, chunksize: int = 500_000) -> pd.DataFrame:
    """EPC 를 청크로 읽어 outcode 별 면적 집계 + 우세 LA 산출."""
    usecols = list(EPC_COLS.values())
    # (outcode, la) -> 면적 리스트를 메모리 절약 위해 부분 집계로 누적
    area_frames = []
    la_counter = {}     # outcode -> {la: count}

    for f in iter_epc_files(epc_path):
        reader = pd.read_csv(
            f, usecols=lambda c: c in usecols, dtype=str,
            chunksize=chunksize, on_bad_lines="skip", low_memory=False,
        )
        for chunk in reader:
            chunk = chunk.rename(columns={v: k for k, v in EPC_COLS.items()})
            chunk["outcode"] = chunk["postcode"].map(to_outcode)
            chunk["floor_area"] = pd.to_numeric(chunk["floor_area"],
                                                errors="coerce")
            chunk = chunk.dropna(subset=["outcode", "floor_area"])
            chunk = chunk[(chunk["floor_area"] >= FLOOR_AREA_MIN)
                          & (chunk["floor_area"] <= FLOOR_AREA_MAX)]
            if chunk.empty:
                continue
            area_frames.append(chunk[["outcode", "floor_area"]].copy())
            # outcode 별 LA 빈도 누적 (우세 LA 결정용)
            for oc, la in zip(chunk["outcode"], chunk["la"].fillna("")):
                if not la:
                    continue
                d = la_counter.setdefault(oc, {})
                d[la] = d.get(la, 0) + 1

    if not area_frames:
        raise SystemExit("EPC 에서 유효한 레코드를 찾지 못했습니다. 경로/스키마 확인.")

    areas = pd.concat(area_frames, ignore_index=True)
    agg = (areas.groupby("outcode")["floor_area"]
                .agg(median_floor_area="median", n_epc="size")
                .reset_index())

    # outcode -> 최빈 LA
    dom_la = {oc: max(d, key=d.get) for oc, d in la_counter.items() if d}
    agg["local_authority"] = agg["outcode"].map(dom_la)
    return agg


# ----------------------------------------------------------------------
# 2) 임대료 로딩 (LA 단위)
# ----------------------------------------------------------------------
def load_rent(rent_csv: str) -> pd.DataFrame:
    """
    정규화된 임대료 CSV 로드.
    필수 컬럼: la_code, median_rent  (선택: bedrooms)
    bedrooms 가 있으면 LA 별 전체 대표값으로 평균(가중X, 단순 median) 후 사용.
    """
    df = pd.read_csv(rent_csv, dtype={"la_code": str})
    df.columns = [c.strip().lower() for c in df.columns]
    if "la_code" not in df.columns or "median_rent" not in df.columns:
        raise SystemExit("rent CSV 에는 최소 la_code, median_rent 컬럼이 필요합니다.")
    df["median_rent"] = pd.to_numeric(df["median_rent"], errors="coerce")
    df = df.dropna(subset=["la_code", "median_rent"])
    # 침실수별이면 LA 대표값(중위)으로 축약
    rent_la = (df.groupby("la_code")["median_rent"].median()
                 .reset_index())
    return rent_la


# ----------------------------------------------------------------------
# 3) outcode 마스터 목록
# ----------------------------------------------------------------------
def load_master_outcodes(onspd_csv: str | None,
                          fallback: pd.Series) -> pd.Index:
    """ONSPD 가 있으면 거기서 모든 outcode 를, 없으면 EPC outcode 를 마스터로."""
    if onspd_csv:
        # ONSPD 는 큰 파일 -> 우편번호 컬럼만 청크로 읽음
        seen = set()
        reader = pd.read_csv(onspd_csv, usecols=lambda c: c.lower() in ("pcds", "pcd"),
                             dtype=str, chunksize=1_000_000, low_memory=False)
        for chunk in reader:
            col = "pcds" if "pcds" in [c.lower() for c in chunk.columns] else "pcd"
            chunk.columns = [c.lower() for c in chunk.columns]
            seen.update(chunk[col].map(to_outcode).dropna().unique())
        return pd.Index(sorted(seen), name="outcode")
    return pd.Index(sorted(fallback.dropna().unique()), name="outcode")


# ----------------------------------------------------------------------
# 4) 옵션: 매물 스크래핑 데이터로 outcode 단위 실제 임대료 augment
# ----------------------------------------------------------------------
def load_listings(listings_csv: str) -> pd.DataFrame:
    """
    스크래핑한 매물 CSV 로 outcode 단위 실거래(호가) 평당가를 직접 계산.
    필수 컬럼: postcode, monthly_rent, floor_area_sqm
    (floor_area 가 sqft 면 미리 ㎡로 변환해 둘 것)
    반환: outcode, rent_per_sqm_listings, n_listings
    """
    df = pd.read_csv(listings_csv)
    df.columns = [c.strip().lower() for c in df.columns]
    df["outcode"] = df["postcode"].map(to_outcode)
    df["monthly_rent"] = pd.to_numeric(df["monthly_rent"], errors="coerce")
    df["floor_area_sqm"] = pd.to_numeric(df["floor_area_sqm"], errors="coerce")
    df = df.dropna(subset=["outcode", "monthly_rent", "floor_area_sqm"])
    df = df[(df["floor_area_sqm"] >= FLOOR_AREA_MIN)
            & (df["floor_area_sqm"] <= FLOOR_AREA_MAX)]
    df["psm"] = df["monthly_rent"] / df["floor_area_sqm"]
    out = (df.groupby("outcode")["psm"]
             .agg(rent_per_sqm_listings="median", n_listings="size")
             .reset_index())
    return out


# ----------------------------------------------------------------------
# 메인 파이프라인
# ----------------------------------------------------------------------
def build(args) -> pd.DataFrame:
    epc = load_epc(args.epc)                       # outcode, median_floor_area, n_epc, local_authority
    rent = load_rent(args.rent)                    # la_code, median_rent

    master = load_master_outcodes(args.onspd, epc["outcode"])
    df = pd.DataFrame(index=master).reset_index()  # 모든 outcode 한 행씩

    df = df.merge(epc, on="outcode", how="left")
    df = df.merge(rent, left_on="local_authority", right_on="la_code", how="left")

    # 핵심 계산: 단위면적당 임대료
    df["rent_per_sqm"] = df["median_rent"] / df["median_floor_area"]

    # 선택: 매물 스크래핑 값으로 덮어쓰기/보강 (있을 때 우선)
    if args.listings:
        lst = load_listings(args.listings)
        df = df.merge(lst, on="outcode", how="left")
        # 스크래핑 표본이 충분하면 그 값을 우선 사용
        use_listings = df["n_listings"].fillna(0) >= MIN_SAMPLE
        df.loc[use_listings, "rent_per_sqm"] = df.loc[use_listings,
                                                      "rent_per_sqm_listings"]

    df["rent_per_sqft"] = df["rent_per_sqm"] * SQM_PER_SQFT
    df["rent_per_pyeong"] = df["rent_per_sqm"] * SQM_PER_PYEONG
    df["low_sample"] = df["n_epc"].fillna(0) < MIN_SAMPLE

    cols = ["outcode", "rent_per_sqm", "rent_per_sqft", "rent_per_pyeong",
            "median_rent", "median_floor_area", "local_authority",
            "n_epc", "low_sample"]
    cols = [c for c in cols if c in df.columns]
    df = df[cols].sort_values("outcode").reset_index(drop=True)
    return df


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="UK outcode 단위 단위면적당 임대료 산출")
    p.add_argument("--epc",
                   help="EPC bulk 디렉터리(certificates.csv 포함) 또는 단일 CSV")
    p.add_argument("--rent",
                   help="정규화된 LA 단위 임대료 CSV (la_code, median_rent[, bedrooms])")
    p.add_argument("--onspd", default=None,
                   help="(선택) ONSPD CSV — 모든 outcode 마스터 목록용")
    p.add_argument("--listings", default=None,
                   help="(선택) 매물 스크래핑 CSV (postcode, monthly_rent, floor_area_sqm)")
    p.add_argument("--out", default="uk_outcode_rent_per_area.csv")
    p.add_argument("--selftest", action="store_true",
                   help="합성 데이터로 로직 검증 후 종료")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.selftest:
        return run_selftest()
    if not args.epc or not args.rent:
        sys.exit("실제 실행에는 --epc 와 --rent 가 필요합니다 (또는 --selftest).")
    df = build(args)
    df.to_csv(args.out, index=False, encoding="utf-8-sig")
    n_total = len(df)
    n_valued = df["rent_per_sqm"].notna().sum()
    print(f"[완료] {args.out}")
    print(f"  - 총 outcode 행 수 : {n_total}")
    print(f"  - 임대료 산출된 행 : {n_valued} ({n_valued/n_total:.1%})")
    print(df.head(8).to_string(index=False))


# ----------------------------------------------------------------------
# 합성 데이터 자체 검증 (실제 다운로드 없이 로직 확인)
# ----------------------------------------------------------------------
def run_selftest():
    import tempfile, shutil
    tmp = tempfile.mkdtemp(prefix="rentpsm_selftest_")
    try:
        # 가짜 EPC: 두 LA(E001=Bath, E002=Bristol), 세 outcode
        epc_dir = os.path.join(tmp, "epc", "Bath")
        os.makedirs(epc_dir)
        epc_rows = []
        rng = np.random.default_rng(0)
        spec = [("BA2", "E001", 80), ("BA1", "E001", 70), ("BS1", "E002", 60)]
        for oc, la, base in spec:
            for _ in range(50):
                fa = float(base + rng.normal(0, 5))
                epc_rows.append({"POSTCODE": f"{oc} 7AY",
                                 "TOTAL_FLOOR_AREA": round(fa, 1),
                                 "LOCAL_AUTHORITY": la,
                                 "NUMBER_HABITABLE_ROOMS": 3,
                                 "PROPERTY_TYPE": "Flat"})
        pd.DataFrame(epc_rows).to_csv(os.path.join(epc_dir, "certificates.csv"),
                                      index=False)
        # 가짜 임대료: LA 단위
        rent_csv = os.path.join(tmp, "rent.csv")
        pd.DataFrame({"la_code": ["E001", "E002"],
                      "median_rent": [1200, 1500]}).to_csv(rent_csv, index=False)
        # 실행
        args = parse_args(["--epc", os.path.join(tmp, "epc"),
                           "--rent", rent_csv,
                           "--out", os.path.join(tmp, "out.csv")])
        df = build(args)
        print("=== SELF-TEST 결과 ===")
        print(df.to_string(index=False))
        # 검증
        assert {"outcode", "rent_per_sqm"}.issubset(df.columns), "필수 컬럼 누락"
        assert set(df["outcode"]) == {"BA1", "BA2", "BS1"}, "outcode 집합 오류"
        ba2 = df.loc[df.outcode == "BA2"].iloc[0]
        expected = 1200 / ba2["median_floor_area"]
        assert abs(ba2["rent_per_sqm"] - expected) < 1e-6, "psm 계산 오류"
        assert df["rent_per_sqm"].notna().all(), "임대료 결측 발생"
        print("\n[OK] 필수 컬럼 존재, outcode 매핑/계산 검증 통과")
        return 0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
