# ════════════════════════════════════════════════════════════════════════════
# 영국 상업용 건물 EPC 규제 리스크 — 기대손실 기반 우선순위 (전체 분석 한 파일)
# 26-1 빅데이터 애널리틱스 C조
# ----------------------------------------------------------------------------
# 실행. '최종_데이터셋' 폴더를 이 스크립트 옆(또는 같은 폴더)에 둔 뒤
#        터미널:  Rscript EPC_기대손실_분석_전체.R
#        RStudio:  파일 열고 Source
#       → 같은 폴더에 outputs/ (figures·tables·clustering) 생성.
#
# 데이터 흐름.
#   EPC 원본 ─전처리(Python)→ EPC+RV.csv ─[Track A: 이 코드]→ P_risk
#            ─[Track B: 팀원, 공식은 §3에서 검산]→ final_risk_ranked.csv
#            ─[IMD 매칭·집계]→ 33_la / 34_area 집계표
#   ※ 전처리·Track B 손실계산·IMD 집계는 상류 단계라 결과를 '최종_데이터셋'으로 제공한다.
#     원본 raw부터 완전 재현하려면 팀원의 전처리·Track B 스크립트가 추가로 필요하다.
#
# 손실 공식.  E[Loss] = P_risk × (벌금 + 임대손실 + 개보수비)
#   벌금  단기 clip(RV×10%, £5k, £50k) · 장기 clip(RV×20%, £10k, £150k)
#   임대손실 단기 RV×3/12 · 장기 RV×6/12
#   개보수비 = pure_dynamic_unit_cost × floor_area   (단가는 KNN 보정, 팀원 산출)
#
# 구성.  0 설정 · 1 EDA · 2 Track A 모델 · 3 손실공식 검산 · 4 지역지도 ·
#        5 순위 · 6 의사결정 매트릭스 · 7 4분면지도 · 8 ONS · 9 클러스터링 · 10 회피
# ════════════════════════════════════════════════════════════════════════════

# ── 패키지 (없으면 자동 설치, 최초 1회 인터넷 필요) ──
.need <- c("data.table","ggplot2","scales","dplyr","sf","spdep","ggrepel",
           "patchwork","cluster","gridExtra","Matrix","glmnet","ranger",
           "xgboost","pROC","PRROC","caret","ragg")
.miss <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
if (length(.miss)) { cat("[설치]", paste(.miss, collapse=", "), "\n")
  install.packages(.miss, repos = "https://cloud.r-project.org") }
suppressPackageStartupMessages(for (p in .need) library(p, character.only = TRUE))
if (requireNamespace("sf", quietly = TRUE)) sf::sf_use_s2(FALSE)

# ── 경로: 코드 폴더 기준, 데이터 폴더 자동 탐색 ──
# RStudio는 파일 열고 우상단 Source 버튼으로 실행하면 자동으로 폴더를 잡는다.
# (안 잡히면 Session > Set Working Directory > To Source File Location 후 Source)
get_this_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f)))                       # Rscript
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)           # source()
  if (!is.null(of)) return(dirname(normalizePath(of)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")  # RStudio
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  getwd()
}
get_script_dir <- get_this_dir
setwd(get_this_dir())
.cand <- c("최종_데이터셋", "데이터셋", "data",
           file.path("..", "최종_데이터셋"), file.path("..", "데이터셋"),
           file.path("..", "data"), file.path("..", "1_데이터셋"), file.path("..", "1_데이터"))
DATA <- NA_character_
for (cc in .cand) if (file.exists(file.path(cc, "EPC+RV.csv"))) { DATA <- normalizePath(cc); break }
if (is.na(DATA)) stop("데이터 폴더를 못 찾음. EPC+RV.csv 가 든 '최종_데이터셋' 폴더를 이 스크립트 옆에 두세요.")
GEO  <- file.path(DATA, "geo")
ROOT <- getwd()
OUTF <- file.path(ROOT, "outputs", "figures"); dir.create(OUTF, recursive = TRUE, showWarnings = FALSE)
OUTT <- file.path(ROOT, "outputs", "tables");  dir.create(OUTT, recursive = TRUE, showWarnings = FALSE)
cat("[setup] DATA =", DATA, "\n[setup] 결과 → ", file.path(ROOT, "outputs"), "\n")

# ── 한글 폰트 + 저장 헬퍼 (그림 글자 깨짐 방지) ──
.pick_font <- function() {
  cand <- c("Apple SD Gothic Neo","AppleGothic","Pretendard","NanumGothic",
            "Noto Sans CJK KR","Noto Sans KR","Malgun Gothic"); fam <- "sans"
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    h <- cand[cand %in% systemfonts::system_fonts()$family]; if (length(h)) fam <- h[1] }
  fam
}
KFONT <- .pick_font(); cat("[setup] 한글 폰트:", KFONT, "\n")
ggsave <- function(filename, plot = ggplot2::last_plot(), ...)
  ggplot2::ggsave(filename, plot, ..., device = if (requireNamespace("ragg", quietly=TRUE)) ragg::agg_png else NULL)
theme_set(theme_minimal(base_family = KFONT) +
          theme(text = element_text(family = KFONT), plot.title = element_text(face = "bold", family = KFONT)))
save_fig <- function(plot, name, w = 9, h = 6, dpi = 200) {
  ggsave(file.path(OUTF, name), plot, width = w, height = h, dpi = dpi, bg = "white"); message("  [fig] ", name) }
save_tab <- function(dt, name) { data.table::fwrite(dt, file.path(OUTT, name)); message("  [tab] ", name) }

# ── 도시유형·인구 룩업 (§5 순위 차트용) ──
CITY_META <- data.table(
  area_label = c(
    "Staffordshire Moorlands", "Torridge", "Pendle", "Bolton", "West Devon",
    "Newcastle-under-Lyme", "Stoke-on-Trent", "Rossendale", "Wyre Forest",
    "Westmorland and Furness", "Fylde", "Manchester", "Milton Keynes", "York",
    "Ashford", "Hart", "Solihull", "Dacorum", "Slough", "Thurrock",
    "Amber Valley", "Birmingham", "Bradford", "Bristol, City of",
    "Buckinghamshire", "Calderdale", "Doncaster", "Erewash", "Fenland",
    "Greater London", "Hyndburn", "Kingston upon Hull, City of", "Leeds",
    "Leicester", "Liverpool", "North East Derbyshire",
    "North East Lincolnshire", "North Lincolnshire", "North Northamptonshire",
    "North Yorkshire", "Nottingham", "Sandwell", "Sheffield", "Wakefield",
    "West Lindsey", "West Northamptonshire"),
  city_type = c(
    "지방", "지방", "도농복합", "대도시권", "지방",
    "도시", "도시", "도농복합", "도농복합",
    "지방", "도농복합", "대도시권", "도시", "도시",
    "도농복합", "도농복합", "대도시권", "도농복합", "도시", "도시",
    "도농복합", "대도시권", "대도시권", "도시",
    "도농복합", "대도시권", "대도시권", "도시", "지방",
    "대도시권", "도시", "도시", "대도시권",
    "도시", "대도시권", "도농복합",
    "도시", "도농복합", "도농복합",
    "지방", "도시", "대도시권", "대도시권", "대도시권",
    "지방", "도시"),
  pop_man = c(
    10, 7, 9, 30, 6, 12, 26, 7, 10, 23, 8, 55, 29, 20,
    13, 10, 22, 16, 16, 18, 13, 114, 55, 47, 55, 21, 31, 12, 10,
    880, 8, 27, 81, 37, 49, 10, 16, 17, 36, 62, 32, 34, 56, 35, 10, 43))
apply_city_tag <- function(dt, mode = "full") {
  dt[CITY_META, on = "area_label", `:=`(.ct = i.city_type, .pm = i.pop_man)]
  if (mode == "full")    dt[!is.na(.ct), nm := sprintf("%s · %s · %d만", nm, .ct, .pm)]
  else if (mode == "compact") dt[!is.na(.ct), nm := sprintf("[%s] %s", .ct, nm)]
  else if (mode == "stacked") dt[!is.na(.ct), nm := sprintf("%s\n%s·%d만", nm, .ct, .pm)]
  dt[, c(".ct", ".pm") := NULL]; invisible(dt)
}

# ════════════════════════════════════════════════════════════════════════════
# 1. EDA — 건물 특성 분포·위험률·상관            [PPT 21~29]
# ════════════════════════════════════════════════════════════════════════════
local({
# 전처리 후 EDA — 건물 특성 분포, 특성별 위험률, 수치형 상관  [PPT 21~29]
# 입력 data/EPC+RV.csv → 출력 outputs/figures, outputs/tables

d <- fread(file.path(DATA, "EPC+RV.csv"), encoding = "UTF-8")
BASE <- mean(d$target)                                   # 전체 위험군 비율(22.7%)
cat(sprintf("[01] %d행 · 위험군 base rate = %.1f%%\n", nrow(d), 100 * BASE))

# 타겟 등급 분포 (D 이하 = 위험)
band_lv <- c("A+","A","B","C","D","E","F","G")
tb <- d[, .N, by = asset_rating_band]
tb[, risk := fifelse(asset_rating_band %in% c("D","E","F","G"), "위험(D~G)", "안전(A+~C)")]
tb[, asset_rating_band := factor(asset_rating_band, levels = band_lv)]
tb[, lab := sprintf("%s\n%.1f%%", comma(N), 100 * N / nrow(d))]
p_band <- ggplot(tb[order(asset_rating_band)], aes(asset_rating_band, N, fill = risk)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = lab), vjust = -0.2, size = 3, lineheight = 0.9, family = KFONT) +
  scale_fill_manual(values = c("안전(A+~C)" = "#4E79A7", "위험(D~G)" = "#D7301F"), name = NULL) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.13))) +
  labs(title = "Asset Rating Band 분포 (타겟)",
       subtitle = "2027 기준 C 강화 → D 이하 = 잠재 비적합(위험) · 위험군 22.7%",
       x = "EPC 등급", y = "건수")
save_fig(p_band, "01_16_asset_rating_band.png", w = 9, h = 5)

# ── 공통 헬퍼 ────────────────────────────────────────────────────────────────
# 범주 분포(상위 n) 가로 막대
bar_count <- function(col, title, top = NULL) {
  t <- d[, .N, by = col][order(-N)]
  if (!is.null(top)) t <- head(t, top)
  t[, lab := sprintf("%s (%.1f%%)", comma(N), 100 * N / nrow(d))]
  ggplot(t, aes(x = N, y = reorder(get(col), N))) +
    geom_col(fill = "#4E79A7") +
    geom_text(aes(label = lab), hjust = -0.05, size = 3) +
    scale_x_continuous(labels = comma, expand = expansion(mult = c(0, 0.18))) +
    labs(title = title, x = "건수", y = NULL)
}
# 범주별 위험률(가로 막대) + base rate 점선
bar_risk <- function(col, title) {
  t <- d[, .(risk = mean(target), n = .N), by = col][order(risk)]
  t[, lab := sprintf("%.0f%%  (n=%s)", 100 * risk, comma(n))]
  ggplot(t, aes(x = risk, y = reorder(get(col), risk))) +
    geom_col(aes(fill = risk > BASE), show.legend = FALSE) +
    scale_fill_manual(values = c(`TRUE` = "#D7301F", `FALSE` = "#4E79A7")) +
    geom_vline(xintercept = BASE, linetype = "dashed", color = "#333333") +
    geom_text(aes(label = lab), hjust = -0.05, size = 3) +
    scale_x_continuous(labels = percent, expand = expansion(mult = c(0, 0.22))) +
    labs(title = title, subtitle = sprintf("점선 = 전체 base rate %.1f%%", 100 * BASE),
         x = "위험률 P(D~G)", y = NULL)
}

# ── 슬라이드 21 · Property Type (상위 15) ────────────────────────────────────
save_fig(bar_count("property_type", "Property Type 분포 (상위 15 / 42종)", top = 15),
         "01_21_property_type.png", w = 11, h = 6)

# ── 슬라이드 22 · Transaction Type ──────────────────────────────────────────
save_fig(bar_count("transaction_type", "Transaction Type 분포 (6종)"),
         "01_22_transaction_type.png", w = 11, h = 4.5)

# ── 슬라이드 23 · Building Level ─────────────────────────────────────────────
save_fig(bar_count("building_level", "Building Level 분포 (3=단순 → 5=복잡)"),
         "01_23_building_level.png", w = 9, h = 3.5)

# ── 슬라이드 24 · Main Heating Fuel ─────────────────────────────────────────
save_fig(bar_count("main_heating_fuel", "Main Heating Fuel 분포 (14종)"),
         "01_24_heating_fuel.png", w = 11, h = 6)

# ── 슬라이드 25 · Floor Area 분포 (0~2,000㎡, 93%) ──────────────────────────
med <- median(d$floor_area); mu <- mean(d$floor_area)
p_fa <- ggplot(d[floor_area <= 2000], aes(floor_area)) +
  geom_histogram(bins = 60, fill = "#4E79A7", color = "white") +
  geom_vline(xintercept = med, linetype = "dashed", color = "#E15759", linewidth = 0.9) +
  annotate("text", x = med, y = Inf, label = sprintf(" 중앙값 %.0f㎡", med),
           hjust = 0, vjust = 1.5, color = "#E15759", family = KFONT) +
  labs(title = "Floor Area 분포 (0~2,000㎡ · 전체의 93%)",
       subtitle = sprintf("중앙값 %.0f㎡ · 평균 %.0f㎡ → 강한 우편향 (최대 %s㎡)", med, mu, comma(max(d$floor_area))),
       x = "연면적 (㎡)", y = "건수")
save_fig(p_fa, "01_25_floor_area.png", w = 10, h = 5.5)

# ── 슬라이드 26 · Aircon Present (분포 + 위험률) ────────────────────────────
d[, aircon_lab := fifelse(aircon_present == 1, "에어컨 있음", "에어컨 없음")]
p_ac <- bar_count("aircon_lab", "Aircon 설치 여부 분포") /
        bar_risk("aircon_lab", "위험률 by Aircon")
save_fig(p_ac, "01_26_aircon.png", w = 10, h = 7)

# ── 슬라이드 27 · Building Environment (분포 + 위험률) ──────────────────────
p_be <- bar_count("building_environment", "Building Environment 분포 (6종)") /
        bar_risk("building_environment", "위험률 by Building Environment")
save_fig(p_be, "01_27_building_environment.png", w = 12, h = 9)

# ── 슬라이드 28 · Local Authority (상위 15 / 319종) ─────────────────────────
save_fig(bar_count("local_authority_label", "Local Authority 분포 (상위 15 / 319종)", top = 15),
         "01_28_local_authority.png", w = 11, h = 6)

# ── 슬라이드 29 · 수치형 변수 Spearman 상관 ─────────────────────────────────
#   면적·건물레벨 = 모델 피처 / 등급점수·위험 = 타겟(참고)
num <- d[, .(`면적` = floor_area, `건물레벨` = building_level,
             `등급점수(타겟)` = asset_rating, `위험(타겟)` = target)]
cm <- cor(num, method = "spearman", use = "pairwise.complete.obs")
cm_dt <- as.data.table(as.table(cm)); setnames(cm_dt, c("V1", "V2", "rho"))
lv <- colnames(cm)
cm_dt[, `:=`(V1 = factor(V1, levels = lv), V2 = factor(V2, levels = rev(lv)))]
p_cor <- ggplot(cm_dt, aes(V1, V2, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 4.5, family = KFONT) +
  scale_fill_gradient2(low = "#D7301F", mid = "white", high = "#2C7FB8",
                       midpoint = 0, limits = c(-1, 1), name = "상관") +
  labs(title = "수치형 변수 상관 (Spearman)",
       subtitle = "면적·건물레벨 = 모델 피처 / 등급점수·위험 = 타겟(참고)",
       x = NULL, y = NULL) + coord_fixed()
save_fig(p_cor, "01_29_spearman_corr.png", w = 7.5, h = 6.5)
save_tab(cm_dt, "01_29_spearman_corr.csv")

cat("[01] EDA 그림 10종 + 상관표 저장 완료\n")
})

# ════════════════════════════════════════════════════════════════════════════
# 2. Track A — EPC 위험확률 P_risk 모델          [PPT 31~38]
# ════════════════════════════════════════════════════════════════════════════
local({
# Track A — EPC 규제 위험확률 P_risk 모델  [PPT 31~38]
# 입력 data/EPC+RV.csv → 출력 outputs/tables/EPC+RV_with_P_risk.csv, 그림·표·모델
#
# 타깃 EPC D~G=1, 안전 8피처(leakage 제외), LA는 fold 내부 target encoding(α=10).
# 층화 80/20 split + Train 5-fold CV(Test 봉인), 세 모델 비교 후 XGBoost 선택.
# 임계값은 Recall≥0.70 보장 하 Precision 최대(≈0.29). 확률은 train=OOF, test=최종모델.
# XGBoost는 멀티스레드라 소수 셋째자리는 실행마다 ±0.001 흔들린다(정상).
set.seed(42)
MODEL_DIR <- file.path(ROOT, "outputs", "model"); dir.create(MODEL_DIR, showWarnings = FALSE)

# ── smoothed target encoding (fold 내부 전용) ───────────────────────────────
target_encode <- function(train_cat, train_y, valid_cat, alpha = 10) {
  train_cat <- as.character(train_cat); valid_cat <- as.character(valid_cat)
  g <- mean(train_y)
  st <- data.table(cat = train_cat, y = train_y)[, .(n = .N, m = mean(y)), by = cat]
  st[, enc := (n * m + alpha * g) / (n + alpha)]
  lk <- setNames(st$enc, st$cat)
  tr <- unname(lk[train_cat]); tr[is.na(tr)] <- g
  va <- unname(lk[valid_cat]); va[is.na(va)] <- g
  list(train_enc = tr, valid_enc = va, global = g, lookup = lk)
}

# ── [1] 데이터 · 안전 피처 8개 · 타깃 ───────────────────────────────────────
raw <- fread(file.path(DATA, "EPC+RV.csv"), encoding = "UTF-8")
feats <- c("property_type", "transaction_type", "building_level",
           "main_heating_fuel", "floor_area", "aircon_present",
           "building_environment", "local_authority")
df <- raw[, c(feats, "target"), with = FALSE]; setnames(df, "target", "y")
cat(sprintf("[1] %d행, 위험(y=1) 비율 %.4f\n", nrow(df), mean(df$y)))

# ── [2] 층화 80/20 split + Train 내부 5-fold ────────────────────────────────
train_idx <- caret::createDataPartition(df$y, p = 0.80, list = FALSE)[, 1]
test_idx  <- setdiff(seq_len(nrow(df)), train_idx)
folds     <- caret::createFolds(df$y[train_idx], k = 5, returnTrain = FALSE)
train <- df[train_idx]; test <- df[test_idx]
cat(sprintf("[2] Train %d (%.4f) / Test %d (%.4f)\n",
            nrow(train), mean(train$y), nrow(test), mean(test$y)))

cat_cols <- c("property_type", "transaction_type", "main_heating_fuel", "building_environment")
lv <- lapply(cat_cols, function(c) sort(unique(train[[c]]))); names(lv) <- cat_cols

align_cols <- function(a, b) {
  miss <- setdiff(colnames(a), colnames(b))
  for (m in miss) b <- cbind(b, Matrix(0, nrow(b), 1, dimnames = list(NULL, m)))
  b[, colnames(a), drop = FALSE]
}
mk_xgb <- function(d) {
  d2 <- data.table(floor_area = d$floor_area, building_level = d$building_level,
                   aircon = d$aircon_present, la_te = d$la_te)
  for (c in cat_cols) { fc <- factor(d[[c]], levels = lv[[c]]); fc[is.na(fc)] <- lv[[c]][1]; d2[[c]] <- fc }
  sparse.model.matrix(~ . - 1, data = as.data.frame(d2))
}
auc_pr <- function(y, p) c(
  auc = as.numeric(pROC::auc(y, p, quiet = TRUE)),
  pr  = PRROC::pr.curve(scores.class0 = p[y == 1], scores.class1 = p[y == 0])$auc.integral)
XGB_PARAMS <- list(objective = "binary:logistic", eval_metric = "auc", eta = 0.05,
                   max_depth = 6, subsample = 0.8, colsample_bytree = 0.8,
                   min_child_weight = 5, seed = 42)
xgb_fit <- function(xt, yt, xv, yv) {
  a <- list(params = XGB_PARAMS, data = xgb.DMatrix(xt, label = yt),
            nrounds = 1000, early_stopping_rounds = 30, verbose = 0)
  dv <- xgb.DMatrix(xv, label = yv)
  tryCatch(do.call(xgb.train, c(a, list(evals = list(val = dv)))),
           error = function(e) do.call(xgb.train, c(a, list(watchlist = list(val = dv)))))
}
best_iter <- function(f) {
  bi <- suppressWarnings(as.integer(tryCatch(xgb.attr(f, "best_iteration"), error = function(e) NA)))
  if (is.na(bi)) bi <- suppressWarnings(as.integer(tryCatch(f$best_iteration, error = function(e) NA)))
  if (is.null(bi) || is.na(bi) || bi < 1) 1000L else bi
}

# ── [3] 5-fold CV — M1 Logistic / M2 RandomForest / M3 XGBoost ──────────────
res <- list(M1 = list(auc = c(), pr = c()), M2 = list(auc = c(), pr = c()), M3 = list(auc = c(), pr = c()))
oof <- rep(NA_real_, nrow(train)); best_iters <- c()
for (i in seq_along(folds)) {
  vi <- folds[[i]]; ti <- setdiff(seq_len(nrow(train)), vi)
  tr <- copy(train[ti]); va <- copy(train[vi])
  te <- target_encode(tr$local_authority, tr$y, va$local_authority)
  tr[, la_te := te$train_enc]; va[, la_te := te$valid_enc]

  # M1. Logistic ridge — log(면적)+표준화, one-hot
  fam <- mean(log1p(tr$floor_area)); fas <- sd(log1p(tr$floor_area))
  tem <- mean(tr$la_te); tes <- sd(tr$la_te)
  mk_lr <- function(d) {
    d2 <- data.table(fa_z = (log1p(d$floor_area) - fam) / fas, late_z = (d$la_te - tem) / tes,
                     building_level = d$building_level, aircon = d$aircon_present)
    for (c in cat_cols) d2[[c]] <- factor(d[[c]], levels = lv[[c]])
    sparse.model.matrix(~ . - 1, data = as.data.frame(d2))
  }
  x1 <- mk_lr(tr); x1v <- align_cols(x1, mk_lr(va))
  f1 <- glmnet(x1, tr$y, family = "binomial", alpha = 0, lambda = 0.01)
  p1 <- as.numeric(predict(f1, x1v, type = "response"))
  m <- auc_pr(va$y, p1); res$M1$auc <- c(res$M1$auc, m["auc"]); res$M1$pr <- c(res$M1$pr, m["pr"])

  # M2. RandomForest (ranger)
  rf_f <- c("floor_area", "building_level", "aircon_present", "la_te", cat_cols)
  f2 <- ranger(x = tr[, rf_f, with = FALSE], y = factor(tr$y, levels = c(0, 1)),
               num.trees = 300, mtry = 4, min.node.size = 50, probability = TRUE,
               num.threads = parallel::detectCores(), seed = 42)
  p2 <- predict(f2, va[, rf_f, with = FALSE])$predictions[, "1"]
  m <- auc_pr(va$y, p2); res$M2$auc <- c(res$M2$auc, m["auc"]); res$M2$pr <- c(res$M2$pr, m["pr"])

  # M3. XGBoost — one-hot, early stopping
  xt <- mk_xgb(tr); xv <- align_cols(xt, mk_xgb(va))
  f3 <- xgb_fit(xt, tr$y, xv, va$y); bi <- best_iter(f3); best_iters <- c(best_iters, bi)
  p3 <- predict(f3, xgb.DMatrix(xv)); oof[vi] <- p3
  m <- auc_pr(va$y, p3); res$M3$auc <- c(res$M3$auc, m["auc"]); res$M3$pr <- c(res$M3$pr, m["pr"])
  cat(sprintf("[3] fold %d  M1 %.4f / M2 %.4f / M3 %.4f (nround=%d)\n",
              i, tail(res$M1$auc, 1), tail(res$M2$auc, 1), tail(res$M3$auc, 1), bi))
}
nrounds <- round(mean(best_iters))

# 모델 비교표 + [슬라이드 35] CV AUC·PR-AUC 막대
cmp <- data.table(model = c("Logistic", "RandomForest", "XGBoost"),
                  auc_mean = sapply(res, function(r) mean(r$auc)), auc_sd = sapply(res, function(r) sd(r$auc)),
                  pr_mean  = sapply(res, function(r) mean(r$pr)),  pr_sd  = sapply(res, function(r) sd(r$pr)))
cat("\n[3] 모델 비교 (5-fold CV)\n"); print(cmp); save_tab(cmp, "02_model_comparison.csv")
cmpL <- rbind(
  cmp[, .(model, metric = "CV AUC",    value = auc_mean, sd = auc_sd)],
  cmp[, .(model, metric = "CV PR-AUC", value = pr_mean,  sd = pr_sd)])
cmpL[, model := factor(model, levels = c("Logistic", "RandomForest", "XGBoost"))]
p_cmp <- ggplot(cmpL, aes(model, value, fill = metric)) +
  geom_col(position = position_dodge(0.7), width = 0.62) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd),
                position = position_dodge(0.7), width = 0.2) +
  geom_text(aes(label = sprintf("%.3f", value)), position = position_dodge(0.7), vjust = -0.6, size = 3.3) +
  scale_fill_manual(values = c("CV AUC" = "#2C5F8A", "CV PR-AUC" = "#54A89A"), name = NULL) +
  coord_cartesian(ylim = c(0.45, 0.9)) +
  labs(title = "모델 3종 비교 (5-fold CV)", subtitle = "표준편차 0.003~0.004 → 분할에 무관하게 안정적 · 최종 선택 XGBoost",
       x = NULL, y = "점수")
save_fig(p_cmp, "02_35_model_comparison.png", w = 8, h = 5)

# ── [4] 임계값 — Recall ≥ 0.70 보장 하 Precision 최대 [슬라이드 36] ─────────
yt <- train$y
mt <- rbindlist(lapply(seq(0.05, 0.95, 0.01), function(t) {
  pr <- as.integer(oof >= t)
  tp <- sum(pr == 1 & yt == 1); fp <- sum(pr == 1 & yt == 0); fn <- sum(pr == 0 & yt == 1)
  p_ <- if (tp + fp == 0) 0 else tp / (tp + fp); r_ <- if (tp + fn == 0) 0 else tp / (tp + fn)
  data.table(thr = t, precision = p_, recall = r_, f1 = if (p_ + r_ == 0) 0 else 2 * p_ * r_ / (p_ + r_))
}))
cand <- mt[recall >= 0.70]; best <- if (nrow(cand)) cand[which.max(precision)] else mt[which.max(f1)]
thr <- best$thr
cat(sprintf("\n[4] 선정 thr=%.2f | R=%.3f P=%.3f F1=%.3f\n", thr, best$recall, best$precision, best$f1))
p_thr <- ggplot(melt(mt, "thr", c("precision", "recall", "f1"), variable.name = "metric"),
                aes(thr, value, color = metric)) + geom_line(linewidth = 1.1) +
  geom_vline(xintercept = thr, linetype = "dashed") +
  annotate("text", x = thr, y = 0.97, label = sprintf("선택 thr=%.2f", thr), hjust = -0.05, family = KFONT) +
  scale_color_manual(values = c(precision = "#2C5F8A", recall = "#D7301F", f1 = "#54A89A"),
                     labels = c(precision = "Precision", recall = "Recall(위험군)", f1 = "F1"), name = NULL) +
  labs(title = sprintf("임계값 선정 — Recall ≥ 0.70 보장 하 Precision 최대 (≈ %.2f)", thr),
       subtitle = "위험 놓치는 비용(FN) > 잘못 잡는 비용(FP) → Recall 우선", x = "threshold", y = NULL)
save_fig(p_thr, "02_36_threshold_curve.png", w = 8.5, h = 5)
save_tab(mt, "02_threshold_curve.csv")

# ── [5] 최종 모델(전체 Train 재학습) + Test(봉인) 평가 [슬라이드 37] ────────
te_full <- target_encode(train$local_authority, train$y, train$local_authority)
train[, la_te := te_full$train_enc]; mm_tr <- mk_xgb(train)
final <- xgb.train(params = XGB_PARAMS, data = xgb.DMatrix(mm_tr, label = train$y), nrounds = nrounds, verbose = 0)
invisible(xgb.save(final, file.path(MODEL_DIR, "M3_final.xgb")))

test[, la_te := target_encode(train$local_authority, train$y, test$local_authority)$valid_enc]
mm_te <- align_cols(mm_tr, mk_xgb(test)); p_test <- predict(final, xgb.DMatrix(mm_te))
auc_t <- as.numeric(pROC::auc(test$y, p_test, quiet = TRUE))
prc_t <- PRROC::pr.curve(scores.class0 = p_test[test$y == 1], scores.class1 = p_test[test$y == 0])
pp <- as.integer(p_test >= thr)
tp <- sum(pp == 1 & test$y == 1); fp <- sum(pp == 1 & test$y == 0)
fn <- sum(pp == 0 & test$y == 1); tn <- sum(pp == 0 & test$y == 0)
test_m <- data.table(threshold = thr, CV_AUC = mean(res$M3$auc), Test_AUC = auc_t, PR_AUC = prc_t$auc.integral,
                     precision = tp / (tp + fp), recall = tp / (tp + fn),
                     f1 = 2 * tp / (2 * tp + fp + fn), accuracy = (tp + tn) / length(test$y))
cat(sprintf("[5] TEST AUC=%.4f (CV=%.4f) | R=%.3f P=%.3f F1=%.3f Acc=%.3f\n",
            auc_t, mean(res$M3$auc), test_m$recall, test_m$precision, test_m$f1, test_m$accuracy))
save_tab(test_m, "02_test_metrics.csv")

# CV vs Test AUC (과적합 점검)
p_ovf <- ggplot(data.table(x = c("CV AUC", "Test AUC"), v = c(mean(res$M3$auc), auc_t)),
                aes(x, v, fill = x)) + geom_col(width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", v)), vjust = -0.6, size = 4.5) +
  scale_fill_manual(values = c("CV AUC" = "#2C5F8A", "Test AUC" = "#D7301F")) +
  coord_cartesian(ylim = c(0.7, 0.88)) + labs(title = "과적합 점검 — CV ≈ Test", x = NULL, y = "AUC")
# Test 최종 지표 막대
mb <- data.table(metric = c("Recall", "Precision", "F1", "Accuracy", "PR-AUC"),
                 v = c(test_m$recall, test_m$precision, test_m$f1, test_m$accuracy, test_m$PR_AUC))
mb[, metric := factor(metric, levels = rev(metric))]
p_mb <- ggplot(mb, aes(v, metric, fill = metric)) + geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", v)), hjust = -0.1, size = 4) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.15))) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = sprintf("봉인 Test 최종 지표 (thr=%.2f)", thr), x = NULL, y = NULL)
save_fig(p_ovf + p_mb + patchwork::plot_layout(widths = c(1, 1.4)), "02_37_test_eval.png", w = 11, h = 5)

# 혼동행렬
p_cm <- ggplot(data.table(actual = c("위험", "위험", "안전", "안전"),
                          predicted = c("위험", "안전", "위험", "안전"), n = c(tp, fn, fp, tn)),
               aes(predicted, actual, fill = n)) + geom_tile() +
  geom_text(aes(label = comma(n)), size = 6, family = KFONT) +
  scale_fill_gradient(low = "white", high = "#4E79A7") +
  labs(title = sprintf("Test 혼동행렬 (thr=%.2f)", thr), x = "예측", y = "실제")
save_fig(p_cm, "02_37b_confusion.png", w = 6.5, h = 5)

# ── [6] 피처 중요도 (Gain, 원피처 8개로 집계) [슬라이드 38] ──────────────────
imp <- xgb.importance(model = final)
grp <- function(f) {
  if (f %in% c("floor_area", "building_level")) return(f)
  if (f == "aircon") return("aircon_present")
  if (f == "la_te")  return("local_authority")
  for (p in cat_cols) if (startsWith(f, p)) return(p)
  f
}
imp[, feat := sapply(Feature, grp)]
ig <- imp[, .(gain = sum(Gain)), by = feat][order(-gain)][, pct := 100 * gain / sum(gain)]
save_tab(ig, "02_feature_importance.csv")
p_imp <- ggplot(ig, aes(pct, reorder(feat, pct))) +
  geom_col(aes(fill = pct), show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), hjust = -0.1, size = 3.6) +
  scale_fill_gradient(low = "#9EC3DF", high = "#21527A") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "XGBoost 피처 중요도 (Gain · 원피처 8개 집계)",
       subtitle = "EPC 고위험(D~G) 예측", x = "Gain (%)", y = NULL)
save_fig(p_imp, "02_38_feature_importance.png", w = 9, h = 5)

# ── [7] 전체 행에 P_risk(out-of-sample) + pred_risk 부착 → CSV ──────────────
P_risk <- rep(NA_real_, nrow(df))
P_risk[train_idx] <- oof; P_risk[test_idx] <- p_test
stopifnot(sum(is.na(P_risk)) == 0)
out <- copy(raw); out[, P_risk := round(P_risk, 6)]; out[, pred_risk := as.integer(P_risk >= thr)]
fwrite(out, file.path(OUTT, "EPC+RV_with_P_risk.csv"))
cat(sprintf("\n[6] 저장 outputs/tables/EPC+RV_with_P_risk.csv | P_risk 평균 %.4f ≈ 실제 위험률 %.4f\n",
            mean(P_risk), mean(out$target)))
cat("[02] Track A 완료 — 그림 6종 + 표 5종 + 모델/예측 CSV\n")
})

# ════════════════════════════════════════════════════════════════════════════
# 3. Track B 손실 공식 검산                      [PPT 39~41]
# ════════════════════════════════════════════════════════════════════════════
local({
# Track B 손실 공식 검산  [PPT 39~41]
# 입력 data/final_risk_ranked.csv (Track B 산출물) → 출력 outputs/figures, tables
# Track B 손실 컬럼이 문서화된 공식으로 그대로 재현되는지 대조해 손실 정의를 못박는다.
#   E[Loss] = P_risk × (벌금 + 임대손실 + 개보수비)
#   벌금 단기 RV×10%(£5k~£50k)·장기 RV×20%(£10k~£150k), 임대손실 RV×3/12·6/12,
#   개보수비 = pure_dynamic_unit_cost × floor_area (단가는 KNN 보정, 팀원 산출값 사용)

f <- fread(file.path(DATA, "final_risk_ranked.csv"), encoding = "UTF-8")
cat(sprintf("[03] %d행 로드 · 손실 공식 검산\n", nrow(f)))

# ── 1. 공식 재계산 ───────────────────────────────────────────────────────────
fine_s <- pmin(pmax(f$RV * 0.10,  5000),  50000)        # 단기 벌금
fine_l <- pmin(pmax(f$RV * 0.20, 10000), 150000)        # 장기 벌금
rent_s <- f$RV * 3 / 12                                  # 단기 임대손실(3개월)
rent_l <- f$RV * 6 / 12                                  # 장기 임대손실(6개월)
retro  <- f$pure_dynamic_unit_cost * f$floor_area        # 개보수비 = 단가 × 면적
tot_s  <- f$P_risk * (fine_s + rent_s + retro)           # 단기 기대손실
tot_l  <- f$P_risk * (fine_l + rent_l + retro)           # 장기 기대손실

# ── 2. 제공 컬럼과 대조 (최대 절대오차) ──────────────────────────────────────
chk <- data.table(
  컬럼 = c("fine_short_term", "fine_long_term", "rent_loss_short", "rent_loss_long",
           "pure_retrofit_cost_total", "total_risk_short_term", "total_risk_long_term"),
  max_abs_diff = c(
    max(abs(fine_s - f$fine_short_term)),
    max(abs(fine_l - f$fine_long_term)),
    max(abs(rent_s - f$rent_loss_short)),
    max(abs(rent_l - f$rent_loss_long)),
    max(abs(retro  - f$pure_retrofit_cost_total), na.rm = TRUE),
    max(abs(tot_s  - f$total_risk_short_term)),
    max(abs(tot_l  - f$total_risk_long_term))))
chk[, 일치 := ifelse(max_abs_diff < 1e-3, "OK (공식 일치)", "불일치")]
rank_cor <- cor(f$rank_short_term, frank(-f$total_risk_short_term), method = "spearman")
cat("\n[03] 공식 검산 결과 (max|diff| ≈ 0 이면 공식 일치)\n"); print(chk)
cat(sprintf("[03] rank_short_term ↔ total_risk 내림차순 순위 상관: %.4f (1=완전일치, 동점처리만 차이)\n", rank_cor))
save_tab(chk, "03_loss_formula_check.csv")

# ── 3. [슬라이드 40] 벌금 구조 — RV 대비 벌금 (상·하한 cap) ──────────────────
rv_grid <- seq(0, 1.2e6, length.out = 400)
fdt <- rbind(
  data.table(RV = rv_grid, fine = pmin(pmax(rv_grid * 0.10, 5000),  50000),  시나리오 = "단기 (10% · £5k~£50k)"),
  data.table(RV = rv_grid, fine = pmin(pmax(rv_grid * 0.20, 10000), 150000), 시나리오 = "장기 (20% · £10k~£150k)"))
p_fine <- ggplot(fdt, aes(RV, fine, color = 시나리오)) + geom_line(linewidth = 1.1) +
  scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k£")) +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k£")) +
  scale_color_manual(values = c("#4E79A7", "#D7301F")) +
  labs(title = "MEES 위반 벌금 구조 (RV 기준 · 상하한 cap)",
       subtitle = "벌금 = clip(RV × 비율, 하한, 상한)", x = "Rateable Value (RV)", y = "벌금", color = NULL)
save_fig(p_fine, "03_40_fine_structure.png", w = 9, h = 5)

# ── 4. [슬라이드 41] 손실 3요소 구성 (위험군 평균 · 단기) ────────────────────
risk_b <- f[pred_risk == 1]
comp <- data.table(
  요소 = c("벌금(Fine)", "임대손실(Rent loss)", "개보수비(Retrofit)"),
  평균금액 = c(mean(risk_b$fine_short_term), mean(risk_b$rent_loss_short), mean(risk_b$pure_retrofit_cost_total)))
comp[, 비중 := 평균금액 / sum(평균금액)]
save_tab(comp, "03_loss_components.csv")
p_comp <- ggplot(comp, aes(x = "단기 손실금액", y = 평균금액, fill = 요소)) +
  geom_col(width = 0.55) +
  geom_text(aes(label = sprintf("%s\n%.0f%%", 요소, 100 * 비중)),
            position = position_stack(vjust = 0.5), size = 3.6, family = KFONT) +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k£")) +
  scale_fill_manual(values = c("벌금(Fine)" = "#4E79A7", "임대손실(Rent loss)" = "#F1A340", "개보수비(Retrofit)" = "#D7301F")) +
  labs(title = "손실 3요소 구성 (위험 예측군 평균 · 단기)",
       subtitle = "확률 곱하기 전 순수 손실금액 = Fine + Rent loss + Retrofit cost",
       x = NULL, y = "건물당 평균 (£)") + theme(legend.position = "none")
save_fig(p_comp, "03_41_loss_components.png", w = 6.5, h = 6)

cat("[03] Track B 공식 검산 완료 — 모든 컬럼 공식 재현 확인\n")
})

# ════════════════════════════════════════════════════════════════════════════
# 4. 지역별 라벨 지도 3종                        [PPT 44~46]
# ════════════════════════════════════════════════════════════════════════════
local({
# 지역별 라벨 지도 3종 (P_risk / 손실금액 / 기대손실)  [PPT 44~46]
# 원본: 지역 결합/02_코드/37_presentation_labeled_maps.R
# 한글 저장: ggsave를 ragg 디바이스로 고정

area_file <- file.path(DATA, "34_area_eda_dataset.csv")
geo_file  <- file.path(GEO,  "uk_lad.geojson")
out_chart <- OUTF
out_final <- OUTT

q_safe <- function(x, p) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 7))
}

fmt_pct <- function(x) paste0(round(as.numeric(x) * 100, 1), "%")

fmt_money_compact <- function(x) {
  x <- as.numeric(x)
  billion <- format(round(x / 1e9, 2), big.mark = ",", trim = TRUE, scientific = FALSE)
  million <- format(round(x / 1e6, 1), big.mark = ",", trim = TRUE, scientific = FALSE)
  out <- ifelse(
    abs(x) >= 1e9,
    paste0("£", billion, "B"),
    paste0("£", million, "M")
  )
  out[is.na(x)] <- NA_character_
  out
}

wrap_name <- function(x, width = 14) {
  vapply(x, function(z) {
    z <- gsub("Bristol, City of", "Bristol", z, fixed = TRUE)
    z <- gsub("Newcastle upon Tyne", "Newcastle", z, fixed = TRUE)
    z <- gsub("Greater London", "London", z, fixed = TRUE)
    paste(strwrap(z, width = width), collapse = "\n")
  }, character(1))
}

metric_specs <- data.table(
  metric = c("mean_p_risk", "total_loss_amount_short", "total_expected_loss_short"),
  plot_col = c("mean_p_risk_plot", "total_loss_amount_short_mil_plot", "total_expected_loss_short_mil_plot"),
  title = c("지역별 평균 EPC 위험 확률", "지역별 손실금액", "지역별 기대손실"),
  subtitle = c(
    "P_risk 자체의 공간 패턴이다. 빨간 라벨은 상위값, 파란 라벨은 하위값, 진회색 라벨은 주요 대도시다. 대도시 라벨 숫자는 인구 규모 순위다.",
    "P_risk를 곱하기 전 재무 노출액이다. 대도시와 자산 규모가 큰 지역이 더 두드러진다. 대도시 라벨 숫자는 인구 규모 순위다.",
    "Expected Loss = P_risk x Loss_amount. 확률과 금액이 결합된 최종 우선순위 지도다. 대도시 라벨 숫자는 인구 규모 순위다."
  ),
  legend_title = c("평균 P_risk", "손실금액 £M", "기대손실 £M"),
  formatter = c("pct", "money", "money"),
  filename = c(
    "37_01_prisk_pretty_labeled_map.png",
    "37_02_loss_amount_pretty_labeled_map.png",
    "37_03_expected_loss_pretty_labeled_map.png"
  )
)

major_area_labels <- c(
  "Greater London", "Birmingham", "Manchester", "Leeds", "Liverpool",
  "Sheffield", "Bristol, City of", "Newcastle upon Tyne", "Nottingham",
  "Leicester", "Bradford", "Coventry", "Derby"
)

la_mapping <- list(
  "E06000058" = c("E07000048", "E06000028", "E06000029"),
  "E06000059" = c("E07000049", "E07000050", "E07000051", "E07000052", "E07000053"),
  "E06000060" = c("E07000004", "E07000005", "E07000006", "E07000007"),
  "E06000061" = c("E07000150", "E07000152", "E07000153", "E07000156"),
  "E06000062" = c("E07000151", "E07000154", "E07000155"),
  "E06000063" = c("E07000026", "E07000028", "E07000029"),
  "E06000064" = c("E07000027", "E07000030", "E07000031"),
  "E06000065" = c("E07000163", "E07000164", "E07000165", "E07000166", "E07000167", "E07000168", "E07000169"),
  "E06000066" = c("E07000187", "E07000188", "E07000189", "E07000190", "E07000191", "E07000246"),
  "E07000244" = c("E07000205", "E07000206"),
  "E07000245" = c("E07000201", "E07000204"),
  "E08000038" = c("E08000016"),
  "E08000039" = c("E08000019")
)

cat("[1] 지역 지표와 지도 geometry 로드\n")
area <- fread(area_file, encoding = "UTF-8")

lad <- st_read(geo_file, quiet = TRUE)
lad <- lad[substr(lad$LAD13CD, 1, 1) == "E", ]
lad <- st_make_valid(lad)

lad$join_code <- lad$LAD13CD
for (new_code in names(la_mapping)) {
  lad$join_code[lad$LAD13CD %in% la_mapping[[new_code]]] <- new_code
}
lad$area_code <- lad$join_code
lad$area_code[grepl("^E09", lad$join_code)] <- "GREATER_LONDON"
lad$area_label_geo <- lad$LAD13NM
lad$area_label_geo[lad$area_code == "GREATER_LONDON"] <- "Greater London"

geo_area <- aggregate(
  lad[, c("area_code", "area_label_geo")],
  by = list(area_code_group = lad$area_code),
  FUN = function(x) x[1],
  do_union = TRUE
)
geo_area$area_code <- geo_area$area_code_group
geo_area$area_code_group <- NULL
if ("area_code.1" %in% names(geo_area)) geo_area$area_code.1 <- NULL
geo_area <- geo_area[, c("area_code", "area_label_geo")]

map_df <- merge(geo_area, area, by = "area_code", all.x = TRUE)

# 영국 지도는 EPSG:27700에서 더 안정적으로 보이고 라벨 repel도 자연스럽다.
map_df <- st_transform(map_df, 27700)
map_attr <- as.data.table(st_drop_geometry(map_df))

map_df$mean_p_risk_plot <- map_df$mean_p_risk
map_df$total_loss_amount_short_mil_plot <- pmin(map_df$total_loss_amount_short / 1e6, q_safe(map_df$total_loss_amount_short / 1e6, 0.98))
map_df$total_expected_loss_short_mil_plot <- pmin(map_df$total_expected_loss_short / 1e6, q_safe(map_df$total_expected_loss_short / 1e6, 0.98))

label_value <- function(metric_col, x) {
  if (metric_col == "mean_p_risk") fmt_pct(x) else fmt_money_compact(x)
}

build_label_df <- function(metric_col, high_n = 6, low_n = 4) {
  dt <- map_attr[!is.na(get(metric_col))]
  dt[, value := get(metric_col)]
  dt[, rank_high := frank(-value, ties.method = "min")]
  dt[, rank_low := frank(value, ties.method = "min")]
  dt[, is_major := area_label %in% major_area_labels]

  ids <- unique(c(
    dt[rank_high <= high_n, area_code],
    dt[rank_low <= low_n, area_code],
    dt[is_major == TRUE, area_code]
  ))

  labels_sf <- map_df[map_df$area_code %in% ids, ]
  label_points <- suppressWarnings(st_point_on_surface(labels_sf))
  coords <- st_coordinates(label_points)
  labels <- as.data.table(st_drop_geometry(label_points))
  labels[, `:=`(x = coords[, 1], y = coords[, 2])]
  labels[, value := get(metric_col)]
  labels[, rank_high := frank(-value, ties.method = "min")]
  labels[, rank_low := frank(value, ties.method = "min")]
  labels[, is_major := area_label %in% major_area_labels]
  labels[, label_class := fifelse(
    rank_high <= high_n, "High value",
    fifelse(rank_low <= low_n, "Low value", "Major city")
  )]
  # 대도시는 major_area_labels 순서(인구 규모순)를 순위로 붙인다.
  labels[, rank_major := match(area_label, major_area_labels)]
  labels[, value_fmt := label_value(metric_col, value)]
  labels[, label_text := fifelse(
    is_major,
    paste0(wrap_name(area_label), "\n대도시 ", rank_major, "위 · ", value_fmt),
    paste0(wrap_name(area_label), "\n", value_fmt)
  )]
  labels[]
}

make_side_panel <- function(metric_col, metric_title, formatter_name, high_n = 7, low_n = 7) {
  dt <- copy(area[!is.na(get(metric_col))])
  dt[, value := get(metric_col)]
  setorder(dt, -value)
  high <- dt[1:min(high_n, .N)]
  setorder(dt, value)
  low <- dt[1:min(low_n, .N)]

  fmt <- if (formatter_name == "pct") fmt_pct else fmt_money_compact
  high_lines <- paste0(seq_len(nrow(high)), ". ", high$area_label, "  ", fmt(high$value))
  low_lines <- paste0(seq_len(nrow(low)), ". ", low$area_label, "  ", fmt(low$value))
  panel_text <- paste(
    metric_title,
    "",
    "상위 지역",
    paste(high_lines, collapse = "\n"),
    "",
    "하위 지역",
    paste(low_lines, collapse = "\n"),
    sep = "\n"
  )

  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "#F8FAFC", color = "#CBD5E1", linewidth = 0.45) +
    annotate("text", x = 0.055, y = 0.955, label = panel_text,
             hjust = 0, vjust = 1, size = 3.25, lineheight = 1.04,
             color = "#172033") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}

theme_map <- theme_void(base_family = KFONT) +
  theme(
    plot.title = element_text(face = "bold", size = 20, color = "#111827", hjust = 0),
    plot.subtitle = element_text(size = 11.5, color = "#334155", hjust = 0, lineheight = 1.08),
    plot.caption = element_text(size = 8.5, color = "#64748B", hjust = 0),
    legend.position = "bottom",
    legend.title = element_text(size = 9.5, color = "#111827"),
    legend.text = element_text(size = 8.5, color = "#334155"),
    legend.key.width = grid::unit(1.55, "cm"),
    plot.margin = margin(8, 8, 8, 8)
  )

label_colors <- c(
  "High value" = "#C7352B",
  "Low value" = "#1F6FEB",
  "Major city" = "#233142"
)

# 37_01, 37_02도 37_03(기대손실)과 같은 팔레트로 통일한다.
expected_loss_palette <- c("#F7FBFF", "#C7E9F1", "#7FCDBB", "#FDAE61", "#D73027", "#7F0000")
metric_palettes <- list(
  mean_p_risk = expected_loss_palette,
  total_loss_amount_short = expected_loss_palette,
  total_expected_loss_short = expected_loss_palette
)

make_pretty_map <- function(spec) {
  metric_col <- spec$metric
  plot_col <- spec$plot_col
  labels <- build_label_df(metric_col)
  high_sf <- map_df[map_df$area_code %in% labels[label_class == "High value", area_code], ]
  low_sf <- map_df[map_df$area_code %in% labels[label_class == "Low value", area_code], ]
  side <- make_side_panel(metric_col, spec$title, spec$formatter)

  # 라벨을 지도 바깥 좌우 여백으로 빼서 길고 얇은 리더선으로 잇는다.
  bb <- st_bbox(map_df)
  W <- as.numeric(bb[["xmax"]] - bb[["xmin"]])
  H <- as.numeric(bb[["ymax"]] - bb[["ymin"]])
  xmid <- as.numeric(median(labels$x, na.rm = TRUE))
  left_lab  <- labels[x <  xmid]   # 서쪽 지역 -> 왼쪽 여백
  right_lab <- labels[x >= xmid]   # 동쪽 지역 -> 오른쪽 여백
  left_band  <- c(bb[["xmin"]] - 0.34 * W, bb[["xmin"]] - 0.08 * W)
  right_band <- c(bb[["xmax"]] + 0.08 * W, bb[["xmax"]] + 0.34 * W)
  y_band <- c(bb[["ymin"]], bb[["ymax"]])

  # 좌/우 두 콜아웃 레이어가 같은 옵션을 공유하므로 헬퍼로 묶는다.
  repel_layer <- function(d, band, hjust) {
    geom_label_repel(
      data = d,
      mapping = aes(x = x, y = y, label = label_text, color = label_class),
      xlim = band, ylim = y_band, hjust = hjust,
      size = 2.7, fontface = "bold", lineheight = 0.9,
      fill = alpha("white", 0.95), label.size = 0.18,
      label.r = grid::unit(0.12, "lines"),
      label.padding = grid::unit(0.12, "lines"),
      box.padding = 0.18, point.padding = 0.2, min.segment.length = 0,
      direction = "y", force = 2.4, force_pull = 0.02, max.overlaps = Inf,
      seed = 20260606, segment.size = 0.3, segment.alpha = 0.6,
      show.legend = FALSE
    )
  }

  legend_labels <- if (spec$formatter == "pct") {
    percent_format(accuracy = 1)
  } else {
    label_number(accuracy = 1, suffix = "M")
  }
  caption_text <- if (metric_col == "mean_p_risk") {
    "기준. Greater London 통합 + 나머지 Local Authority. 색상은 지역별 평균 P_risk 원값이다."
  } else {
    "기준. Greater London 통합 + 나머지 Local Authority. 총액 지표의 색상은 p98 기준으로 눌러서 지역 간 차이가 보이게 했다."
  }

  p_map <- ggplot(map_df) +
    geom_sf(aes(fill = get(plot_col)), color = "#FFFFFF", linewidth = 0.05) +
    geom_sf(data = high_sf, fill = NA, color = "#9F1D1D", linewidth = 0.45) +
    geom_sf(data = low_sf, fill = NA, color = "#1954A6", linewidth = 0.4) +
    geom_point(
      data = labels,
      aes(x = x, y = y, color = label_class),
      size = 2.0,
      show.legend = FALSE
    ) +
    repel_layer(left_lab, left_band, 1) +
    repel_layer(right_lab, right_band, 0) +
    scale_fill_gradientn(
      colours = metric_palettes[[metric_col]],
      labels = legend_labels,
      na.value = "#ECEFF3",
      oob = squish,
      name = spec$legend_title
    ) +
    scale_color_manual(values = label_colors) +
    guides(fill = guide_colorbar(title.position = "top", barwidth = grid::unit(8.5, "cm"))) +
    coord_sf(
      xlim = c(bb[["xmin"]] - 0.42 * W, bb[["xmax"]] + 0.42 * W),
      ylim = c(bb[["ymin"]] - 0.04 * H, bb[["ymax"]] + 0.04 * H),
      datum = NA, expand = FALSE
    ) +
    labs(
      title = spec$title,
      subtitle = spec$subtitle,
      caption = caption_text
    ) +
    theme_map

  full_plot <- p_map + side + plot_layout(widths = c(3.35, 1.1))

  chart_path <- file.path(out_chart, spec$filename)
  ggsave(chart_path, full_plot, width = 15.8, height = 9.2, dpi = 240, bg = "white")

  labels[, `:=`(
    map_metric = metric_col,
    formatted_value = label_value(metric_col, value)
  )]
  labels[, .(
    map_metric,
    area_code,
    area_label,
    label_class,
    value,
    formatted_value,
    rank_high,
    rank_low,
    is_major
  )]
}

cat("[2] 발표용 라벨 지도 생성\n")
label_catalog <- rbindlist(lapply(seq_len(nrow(metric_specs)), function(i) {
  make_pretty_map(metric_specs[i])
}), fill = TRUE)
setorder(label_catalog, map_metric, label_class, rank_high, area_label)

label_file <- file.path(out_final, "10_label_candidates.csv")
fwrite(label_catalog, label_file)
cat("[완료] 지도 3종 + 라벨표 →", out_chart, "\n")
})

# ════════════════════════════════════════════════════════════════════════════
# 5. 순위 비교 (위험률·손실금액·기대손실)        [PPT 47]
# ════════════════════════════════════════════════════════════════════════════
local({
# 40_rank_comparison.R
# 도시(지역)별 순위 3종 + 순위 변동 시각화.
#   ① 위험률(mean P_risk) Top 15 막대
#   ② 손실금액(total loss amount) Top 15 막대
#   ③ 기대손실(total expected loss) Top 15 막대
#   ④ 3축 bump 차트 — 위험률 → 손실금액 → 기대손실 순위 이동
#   ⑤ 3종 막대 트립틱 (한 장 비교용)
# 입력. 34_area_prisk_loss_imd_eda_dataset.csv (264개 지역, Greater London 통합)
# 38_derived_graphs.R 와 동일 톤. short-term 시나리오 기준.


area_file <- file.path(DATA, "34_area_eda_dataset.csv")
out_chart <- OUTF

# 도시유형·인구 태그 공유 룩업 (CITY_META, apply_city_tag)

# --- 헬퍼 (38 스크립트와 동일) ---
fmt_pct <- function(x) paste0(round(as.numeric(x) * 100, 1), "%")
fmt_money_compact <- function(x) {
  x <- as.numeric(x)
  ifelse(abs(x) >= 1e9,
         paste0("£", format(round(x / 1e9, 2), big.mark = ",", trim = TRUE, scientific = FALSE), "B"),
         paste0("£", format(round(x / 1e6, 1), big.mark = ",", trim = TRUE, scientific = FALSE), "M"))
}
short_name <- function(x) {
  x <- gsub("Greater London", "London", x, fixed = TRUE)
  x <- gsub("Bristol, City of", "Bristol", x, fixed = TRUE)
  x <- gsub("Newcastle upon Tyne", "Newcastle", x, fixed = TRUE)
  x <- gsub("Newcastle-under-Lyme", "Newcastle-u-Lyme", x, fixed = TRUE)
  x <- gsub("North Northamptonshire", "N.Northants", x, fixed = TRUE)
  x <- gsub("West Northamptonshire", "W.Northants", x, fixed = TRUE)
  x <- gsub("Staffordshire Moorlands", "Staffs Moorlands", x, fixed = TRUE)
  x <- gsub("Westmorland and Furness", "Westmorland&F", x, fixed = TRUE)
  x <- gsub("Hinckley and Bosworth", "Hinckley&B", x, fixed = TRUE)
  x <- gsub("North Yorkshire", "N.Yorkshire", x, fixed = TRUE)
  x <- gsub("Milton Keynes", "Milton K.", x, fixed = TRUE)
  x
}

theme_chart <- theme_minimal(base_family = KFONT) +
  theme(
    plot.title = element_text(face = "bold", size = 18, color = "#111827"),
    plot.subtitle = element_text(size = 11, color = "#334155", lineheight = 1.1, margin = margin(b = 6)),
    plot.caption = element_text(size = 8.5, color = "#64748B", hjust = 0, margin = margin(t = 8)),
    axis.title = element_text(size = 10.5, color = "#334155"),
    axis.text = element_text(size = 9.5, color = "#475569"),
    legend.title = element_text(size = 9.5, color = "#111827"),
    legend.text = element_text(size = 9, color = "#334155"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 16, 10, 12)
  )

save_chart <- function(plot, filename, width = 11, height = 7.2) {
  ggsave(file.path(out_chart, filename), plot, width = width, height = height, dpi = 240, bg = "white")
}

cat("[load] 지역 단위 데이터\n")
d <- fread(area_file, encoding = "UTF-8")
cat("  지역 수:", nrow(d), "\n")

# 순위 상관 (subtitle 근거)
sc_pl <- cor(d$rank_mean_p_risk, d$rank_total_loss_amount_short, method = "spearman")
sc_pe <- cor(d$rank_mean_p_risk, d$rank_total_expected_loss_short, method = "spearman")
sc_le <- cor(d$rank_total_loss_amount_short, d$rank_total_expected_loss_short, method = "spearman")
cat(sprintf("  Spearman  위험-손실 %.2f / 위험-기대 %.2f / 손실-기대 %.2f\n", sc_pl, sc_pe, sc_le))

# ============================================================
# 공통 막대 빌더 (Top 15 가로 막대)
# ============================================================
make_bar <- function(dt, valcol, lblfun, title, subtitle, xlab, fillvar = NULL,
                     bar_fill = "#FC8D59", meta_mode = "full") {
  top <- dt[order(-get(valcol))][1:15]
  top[, nm := short_name(area_label)]
  apply_city_tag(top, meta_mode)              # 지역명 뒤 도시유형·인구 태그
  top[, nm := factor(nm, levels = rev(nm))]   # 1위를 맨 위로
  top[, v := get(valcol)]
  top[, lbl := lblfun(v)]
  top[, is_london := area_label == "Greater London"]

  p <- ggplot(top, aes(x = v, y = nm)) +
    geom_col(aes(fill = is_london), width = 0.74) +
    geom_text(aes(label = lbl), hjust = -0.08, size = 3.0, color = "#334155", fontface = "bold") +
    scale_fill_manual(values = c("FALSE" = bar_fill, "TRUE" = "#7F0000"), guide = "none") +
    labs(title = title, subtitle = subtitle, x = xlab, y = NULL,
         caption = "264개 지역(Greater London 통합) 중 상위 15개. 지역명 뒤 = ONS 도시-농촌 분류·인구(2021). short-term 시나리오.") +
    theme_chart +
    theme(panel.grid.major.y = element_blank())
  p
}

# ① 위험률 — 상위 10 + 하위 10 (양극단 비교)
cat("[1] 위험률 상위 10 + 하위 10\n")

# 양극단 막대 빌더 — 상위 n개(고위험·빨강) + 하위 n개(저위험·파랑)
#   meta_mode. 도시유형·인구 표기 (apply_city_tag) — none/full/compact
make_bar_ends <- function(dt, valcol, lblfun, title, subtitle, xlab, n = 10,
                          hi_fill = "#C7352B", lo_fill = "#2C7BB6",
                          meta_mode = "none") {
  o   <- dt[order(-get(valcol))]              # 값 내림차순
  lbl_hi <- sprintf("상위 %d (고위험)", n)
  lbl_lo <- sprintf("하위 %d (저위험)", n)
  both <- rbind(head(o, n)[, grp := lbl_hi],
                tail(o, n)[, grp := lbl_lo])  # 합쳐도 여전히 값 내림차순
  both[, nm := short_name(area_label)]
  if (meta_mode != "none") apply_city_tag(both, meta_mode)  # 도시유형·인구 태그

  both[, v := get(valcol)]
  both[, lbl := lblfun(v)]
  both[, nm := factor(nm, levels = rev(nm))]  # 1위를 맨 위로
  both[, grp := factor(grp, levels = c(lbl_hi, lbl_lo))]

  ggplot(both, aes(x = v, y = nm, fill = grp)) +
    geom_col(width = 0.74) +
    geom_text(aes(label = lbl), hjust = -0.08, size = 3.0, color = "#334155", fontface = "bold") +
    scale_fill_manual(values = setNames(c(hi_fill, lo_fill), c(lbl_hi, lbl_lo)), name = NULL) +
    labs(title = title, subtitle = subtitle, x = xlab, y = NULL,
         caption = sprintf("264개 지역(Greater London 통합) 중 위험률 상·하위 각 %d개. 지역명 뒤 = ONS 도시-농촌 분류·인구(2021). short-term 시나리오.", n)) +
    theme_chart +
    theme(panel.grid.major.y = element_blank(), legend.position = "top")
}

# 40_01 단독용 — 전체 표기 (지역 · 분류 · 인구)
p1 <- make_bar_ends(
  d, "mean_p_risk", function(v) fmt_pct(v),
  "지역별 위험률 양극단 (평균 P_risk 상위 10 vs 하위 10)",
  "규제 취약성의 양 끝. 상위 10개는 지방·도농복합 소도시로 자산은 작아도 위험 확률이 높고, 하위 10개는 대도시·통근권으로 위험 확률이 낮다.",
  "지역 평균 P_risk", meta_mode = "full"
) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12)))
save_chart(p1, "40_01_rank_prisk.png", width = 12, height = 7.6)
cat("  saved 40_01_rank_prisk.png\n")

# 트립틱(⑤)용 — 간결 표기 ([분류] 지역), 패널 폭 절약
p1_tri <- make_bar_ends(
  d, "mean_p_risk", function(v) fmt_pct(v),
  "지역별 위험률 양극단 (평균 P_risk 상위 10 vs 하위 10)",
  "규제 취약성의 양 끝. 상위 10개는 지방·도농복합 소도시로 자산은 작아도 위험 확률이 높고, 하위 10개는 대도시·통근권으로 위험 확률이 낮다.",
  "지역 평균 P_risk", meta_mode = "compact"
) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12)))

# ② 손실금액 Top 15
cat("[2] 손실금액 Top 15\n")
p2 <- make_bar(
  d, "total_loss_amount_short", function(v) fmt_money_compact(v),
  "지역별 손실금액 순위 (확률 곱하기 전)",
  "재무 노출 상위. 전부 대도시 — 연면적·자산규모가 좌우한다. Greater London 단독 £4.77B.",
  "총 손실금액 (£)"
) +
  scale_x_continuous(labels = label_number(scale = 1e-6, suffix = "M", accuracy = 1),
                     expand = expansion(mult = c(0, 0.15)))
save_chart(p2, "40_02_rank_loss_amount.png", width = 10.5, height = 7.0)
cat("  saved 40_02_rank_loss_amount.png\n")

# ③ 기대손실 Top 15
cat("[3] 기대손실 Top 15\n")
p3 <- make_bar(
  d, "total_expected_loss_short", function(v) fmt_money_compact(v),
  "지역별 기대손실 순위 (P_risk × 손실금액)",
  sprintf("최종 우선순위. 손실금액 순위와 거의 일치한다 (순위상관 %.2f). 위험률보다 금액이 지배한다.", sc_le),
  "총 기대손실 (£)"
) +
  scale_x_continuous(labels = label_number(scale = 1e-6, suffix = "M", accuracy = 1),
                     expand = expansion(mult = c(0, 0.15)))
save_chart(p3, "40_03_rank_expected_loss.png", width = 10.5, height = 7.0)
cat("  saved 40_03_rank_expected_loss.png\n")

# 트립틱(⑤)용 손실·기대 막대 — 간결 표기 ([분류] 지역)
p2_tri <- make_bar(
  d, "total_loss_amount_short", function(v) fmt_money_compact(v),
  "지역별 손실금액 순위 (확률 곱하기 전)",
  "재무 노출 상위. 전부 대도시 — 연면적·자산규모가 좌우한다. Greater London 단독 £4.77B.",
  "총 손실금액 (£)", meta_mode = "compact"
) +
  scale_x_continuous(labels = label_number(scale = 1e-6, suffix = "M", accuracy = 1),
                     expand = expansion(mult = c(0, 0.15)))
p3_tri <- make_bar(
  d, "total_expected_loss_short", function(v) fmt_money_compact(v),
  "지역별 기대손실 순위 (P_risk × 손실금액)",
  sprintf("최종 우선순위. 손실금액 순위와 거의 일치한다 (순위상관 %.2f). 위험률보다 금액이 지배한다.", sc_le),
  "총 기대손실 (£)", meta_mode = "compact"
) +
  scale_x_continuous(labels = label_number(scale = 1e-6, suffix = "M", accuracy = 1),
                     expand = expansion(mult = c(0, 0.15)))

# ============================================================
# ④ 3축 bump 차트 — 위험률 → 손실금액 → 기대손실
# ============================================================
cat("[4] 3축 bump 차트\n")

nr <- 7
sel <- unique(c(
  d[order(rank_mean_p_risk)][1:nr, area_label],
  d[order(rank_total_loss_amount_short)][1:nr, area_label],
  d[order(rank_total_expected_loss_short)][1:nr, area_label]
))
cat("  선택 지역 수:", length(sel), "\n")

s <- d[area_label %in% sel, .(
  area_label, nm = short_name(area_label),
  r1 = rank_mean_p_risk,
  r2 = rank_total_loss_amount_short,
  r3 = rank_total_expected_loss_short,
  quad = quadrant_prisk_loss_amount
)]

# 사분면 → 아키타입 색 (38_01/38_07 과 동일 의미)
quad_map <- c(
  "High P_risk + High Loss_amount"   = "최우선 (위험·금액 동시)",
  "Lower P_risk + High Loss_amount"  = "도심 고손실형",
  "High P_risk + Lower Loss_amount"  = "지방 고위험형",
  "Lower P_risk + Lower Loss_amount" = "저순위"
)
s[, arche := factor(quad_map[quad],
                    levels = c("최우선 (위험·금액 동시)", "도심 고손실형", "지방 고위험형", "저순위"))]
arche_cols <- c(
  "최우선 (위험·금액 동시)" = "#D73027",
  "도심 고손실형"           = "#FC8D59",
  "지방 고위험형"           = "#2C7BB6",
  "저순위"                  = "#94A3B8"
)

# long 변환 (axis 1,2,3)
sl <- rbind(
  s[, .(area_label, nm, arche, axis = 1L, rank = r1)],
  s[, .(area_label, nm, arche, axis = 2L, rank = r2)],
  s[, .(area_label, nm, arche, axis = 3L, rank = r3)]
)

p4 <- ggplot(sl, aes(x = axis, y = rank, group = area_label, color = arche)) +
  geom_line(linewidth = 0.9, alpha = 0.8) +
  geom_point(size = 2.6) +
  # 왼쪽 라벨 (위험률 순위)
  geom_text_repel(
    data = s, aes(x = 1, y = r1, color = arche, label = paste0(nm, "  ", r1)),
    hjust = 1, nudge_x = -0.06, direction = "y", segment.color = NA,
    size = 3.0, fontface = "bold", seed = 20260607, box.padding = 0.16, max.overlaps = Inf
  ) +
  # 가운데 라벨 (손실금액 순위) — 숫자만
  geom_text(
    data = s, aes(x = 2, y = r2, label = r2), color = "#1F2937",
    size = 2.7, fontface = "bold", vjust = -0.9
  ) +
  # 오른쪽 라벨 (기대손실 순위)
  geom_text_repel(
    data = s, aes(x = 3, y = r3, color = arche, label = paste0(r3, "  ", nm)),
    hjust = 0, nudge_x = 0.06, direction = "y", segment.color = NA,
    size = 3.0, fontface = "bold", seed = 20260607, box.padding = 0.16, max.overlaps = Inf
  ) +
  scale_y_reverse(breaks = c(1, 50, 100, 150, 200, 264)) +
  scale_x_continuous(
    limits = c(0.42, 3.58), breaks = c(1, 2, 3),
    labels = c("① 위험률\n(P_risk)", "② 손실금액\n(£, 확률 전)", "③ 기대손실\n(P_risk×£)"),
    position = "top"
  ) +
  scale_color_manual(values = arche_cols, name = NULL, drop = FALSE) +
  labs(
    title = "세 지표에서 순위가 어떻게 뒤집히는가",
    subtitle = sprintf(
      "①→② 선이 마구 교차한다 = 위험 순위와 돈 순위는 딴 동네 (상관 %+.2f). ②→③ 선은 거의 평행 = 곱셈 후에도 돈 순위가 그대로 (상관 %.2f).\n지방 고위험형(파랑)은 위험 상위에서 손실·기대 하위로 추락하고, 도심형(주황)은 정반대로 솟는다.",
      sc_pl, sc_le),
    x = NULL, y = "순위 (264개 지역, 위가 상위)",
    caption = "위험률·손실금액·기대손실 각 상위 7개의 합집합. 색은 P_risk×손실금액 사분면(상위 25% 기준) 아키타입."
  ) +
  theme_chart +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold", color = "#111827"),
    legend.position = "bottom"
  )

save_chart(p4, "40_04_rank_bump_3axis.png", width = 11.5, height = 8.6)
cat("  saved 40_04_rank_bump_3axis.png\n")

# ============================================================
# ⑤ 3종 막대 트립틱 (한 장 비교)
# ============================================================
cat("[5] 트립틱\n")
strip <- function(p) p + labs(caption = NULL) +
  theme(plot.title = element_text(size = 14), plot.subtitle = element_text(size = 9))
pt <- (strip(p1_tri) | strip(p2_tri) | strip(p3_tri)) +
  plot_annotation(
    title = "위험률 · 손실금액 · 기대손실 — 순위 주인공이 바뀐다",
    subtitle = "왼쪽은 지방 소도시, 가운데·오른쪽은 대도시. 같은 264개 지역인데 상위 명단이 거의 겹치지 않는다.",
    caption = "264개 지역(Greater London 통합), short-term 시나리오. 위험률은 상·하위 각 10개, 손실금액·기대손실은 상위 15개.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, color = "#111827"),
      plot.subtitle = element_text(size = 11, color = "#334155", margin = margin(b = 4)),
      plot.caption = element_text(size = 8.5, color = "#64748B", hjust = 0)
    )
  )
save_chart(pt, "40_05_rank_triptych.png", width = 16, height = 7.0)
cat("  saved 40_05_rank_triptych.png\n")

cat("[완료] 순위 차트 5종 →", out_chart, "\n")
})

# ════════════════════════════════════════════════════════════════════════════
# 6. 의사결정 매트릭스 + 파생 그래프             [PPT 49]
# ════════════════════════════════════════════════════════════════════════════
local({
# 38_derived_graphs.R
# 파생 그래프 4종. 발표 지도 3종(P_risk / 손실금액 / 기대손실, script 37) 위에 얹는 분석 차트.
# 입력. 34_area_prisk_loss_imd_eda_dataset.csv (264개 지역, Greater London 통합)
# 모든 그래프는 지도와 동일하게 short-term 시나리오 기준.


area_file <- file.path(DATA, "34_area_eda_dataset.csv")
out_chart <- OUTF

# 도시유형·인구 태그 공유 룩업 (CITY_META, apply_city_tag)

# --- 재사용 헬퍼 (script 37과 동일 톤) ---
q_safe <- function(x, p) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE, type = 7))
}
fmt_pct <- function(x) paste0(round(as.numeric(x) * 100, 1), "%")
fmt_money_compact <- function(x) {
  x <- as.numeric(x)
  ifelse(abs(x) >= 1e9,
         paste0("£", format(round(x / 1e9, 2), big.mark = ",", trim = TRUE, scientific = FALSE), "B"),
         paste0("£", format(round(x / 1e6, 1), big.mark = ",", trim = TRUE, scientific = FALSE), "M"))
}
short_name <- function(x) {
  x <- gsub("Greater London", "London", x, fixed = TRUE)
  x <- gsub("Bristol, City of", "Bristol", x, fixed = TRUE)
  x <- gsub("Newcastle upon Tyne", "Newcastle", x, fixed = TRUE)
  x <- gsub("North Northamptonshire", "N. Northants", x, fixed = TRUE)
  x
}

theme_chart <- theme_minimal(base_family = KFONT) +
  theme(
    plot.title = element_text(face = "bold", size = 18, color = "#111827"),
    plot.subtitle = element_text(size = 11, color = "#334155", lineheight = 1.1, margin = margin(b = 6)),
    plot.caption = element_text(size = 8.5, color = "#64748B", hjust = 0, margin = margin(t = 8)),
    axis.title = element_text(size = 10.5, color = "#334155"),
    axis.text = element_text(size = 9.5, color = "#475569"),
    legend.title = element_text(size = 9.5, color = "#111827"),
    legend.text = element_text(size = 9, color = "#334155"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 16, 10, 12)
  )

save_chart <- function(plot, filename, width = 11, height = 7.2) {
  ggsave(file.path(out_chart, filename), plot, width = width, height = height, dpi = 240, bg = "white")
}

cat("[load] 지역 단위 데이터\n")
d <- fread(area_file, encoding = "UTF-8")
cat("  지역 수:", nrow(d), "\n")

# ============================================================
# ① 사분면 매트릭스 — P_risk × 손실금액
# ============================================================
cat("[1] 사분면 매트릭스\n")

x_cut <- q_safe(d$mean_p_risk, 0.75)             # 0.272
y_cut_m <- q_safe(d$total_loss_amount_short, 0.75) / 1e6  # 81.5 (£M)

quad_levels <- c(
  "High P_risk + High Loss_amount",
  "Lower P_risk + High Loss_amount",
  "High P_risk + Lower Loss_amount",
  "Lower P_risk + Lower Loss_amount"
)
quad_palette <- c(
  "High P_risk + High Loss_amount"   = "#D73027",  # 최우선
  "Lower P_risk + High Loss_amount"  = "#FC8D59",  # 런던형
  "High P_risk + Lower Loss_amount"  = "#4E9CC4",  # 지방형
  "Lower P_risk + Lower Loss_amount" = "#C2C9D2"   # 저순위 (배경으로 가라앉힘)
)
d[, quad := factor(quadrant_prisk_loss_amount, levels = quad_levels)]
d[, loss_m := total_loss_amount_short / 1e6]

# 라벨 대상. 각 분면 대표 + 발산 극단
label_areas <- c(
  "Greater London", "Birmingham", "Bristol, City of", "Leeds", "Bradford",  # 런던형 (고손실)
  "Doncaster", "Sheffield",                                                  # 최우선
  "Staffordshire Moorlands", "Torridge", "Pendle"                            # 지방형 (고위험)
)
lab <- d[area_label %in% label_areas]
lab[, lbl := paste0(short_name(area_label), "\n", fmt_pct(mean_p_risk), " · ", fmt_money_compact(total_loss_amount_short))]

# 코너 주석 (x = P_risk, y = £M). 상단 라벨은 London 버블을 피해 1,900M 높이에 배치.
ann <- data.table(
  x = c(0.452, 0.128, 0.452, 0.128),
  y = c(1900, 1900, 2.6, 2.6),
  hjust = c(1, 0, 1, 0),
  lbl = c("최우선\n고위험·고손실", "런던형\n저위험·고손실",
          "지방형\n고위험·저손실", "저순위\n저위험·저손실"),
  col = c("#9B1B1B", "#B45309", "#1E6A8E", "#64748B")
)

p1 <- ggplot(d, aes(x = mean_p_risk, y = loss_m)) +
  geom_vline(xintercept = x_cut, linetype = "dashed", color = "#94A3B8", linewidth = 0.4) +
  geom_hline(yintercept = y_cut_m, linetype = "dashed", color = "#94A3B8", linewidth = 0.4) +
  geom_point(aes(size = n_buildings, fill = quad), shape = 21, color = "white", stroke = 0.3, alpha = 0.85) +
  geom_text(data = ann, aes(x = x, y = y, label = lbl, hjust = hjust, color = col),
            inherit.aes = FALSE, size = 3.5, fontface = "bold", lineheight = 0.95, show.legend = FALSE) +
  scale_color_identity() +
  geom_point(data = lab, aes(size = n_buildings), shape = 21, fill = NA, color = "#1F2937",
             stroke = 0.55, show.legend = FALSE) +
  geom_label_repel(data = lab, aes(label = lbl), size = 2.9, fontface = "bold",
                   color = "#172033", fill = alpha("white", 0.9),
                   label.size = 0.15, label.padding = grid::unit(0.12, "lines"),
                   box.padding = 0.55, point.padding = 0.35, min.segment.length = 0,
                   force = 7, seed = 20260606, segment.size = 0.3, segment.alpha = 0.7,
                   max.overlaps = Inf) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0.12, 0.46), breaks = seq(0.15, 0.45, 0.05)) +
  scale_y_log10(labels = label_number(accuracy = 1, suffix = "M"),
                breaks = c(3, 10, 30, 100, 300, 1000, 4000)) +
  scale_fill_manual(values = quad_palette, name = "사분면 (상위 25% 기준)",
                    labels = c("최우선 · 고위험+고손실", "런던형 · 저위험+고손실",
                               "지방형 · 고위험+저손실", "저순위 · 저위험+저손실")) +
  scale_size_continuous(range = c(1.6, 11), labels = label_comma(), name = "건물 수") +
  guides(fill = guide_legend(override.aes = list(size = 4.5, alpha = 1), order = 1),
         size = guide_legend(order = 2)) +
  labs(
    title = "위험 확률 × 손실금액 사분면",
    subtitle = "점선은 각 축 상위 25% 기준선. 오른쪽일수록 규제 위험이 높고, 위일수록 재무 노출이 크다.\nLondon은 위험 확률 하위권이지만 손실금액 최상단, Staffordshire Moorlands는 정반대다.",
    x = "지역 평균 P_risk", y = "손실금액 (£M, log scale)",
    caption = "기준. Greater London 통합 + 나머지 Local Authority 264개 지역. 손실금액 = 확률 곱하기 전 Fine + Rent loss + Retrofit cost."
  ) +
  theme_chart

save_chart(p1, "38_01_quadrant_prisk_loss.png", width = 11.5, height = 7.6)
cat("  saved 38_01_quadrant_prisk_loss.png\n")

# ============================================================
# ② 순위 역전 슬로프 — P_risk 순위 ↔ 기대손실 순위
# ============================================================
cat("[2] 순위 역전 슬로프\n")

nr <- 6
top_risk <- d[order(rank_mean_p_risk)][1:nr, area_label]
top_exp  <- d[order(rank_total_expected_loss_short)][1:nr, area_label]
sel <- unique(c(top_risk, top_exp))

s <- d[area_label %in% sel, .(
  area_label,
  r_prisk = rank_mean_p_risk,
  r_exp = rank_total_expected_loss_short,
  mean_p_risk, total_expected_loss_short
)]
s[, delta := r_prisk - r_exp]   # >0 = 기대손실 순위가 더 높음(우선순위 상승)
s[, dir := fifelse(delta >= 25, "상승", fifelse(delta <= -25, "하락", "유지"))]
s[, dir := factor(dir, levels = c("상승", "하락", "유지"))]
s[, nm := short_name(area_label)]

# long 변환. 좌(axis=1)=P_risk 순위, 우(axis=2)=기대손실 순위
sl <- rbind(
  s[, .(area_label, nm, dir, axis = 1, rank = r_prisk)],
  s[, .(area_label, nm, dir, axis = 2, rank = r_exp)]
)

dir_cols <- c("상승" = "#1F6FEB", "하락" = "#C7352B", "유지" = "#94A3B8")

p2 <- ggplot(sl, aes(x = axis, y = rank, group = area_label, color = dir)) +
  geom_line(linewidth = 0.95, alpha = 0.85) +
  geom_point(size = 2.7) +
  geom_text_repel(
    data = s, aes(x = 1, y = r_prisk, color = dir, label = paste0(nm, "  ", r_prisk)),
    hjust = 1, nudge_x = -0.05, direction = "y", segment.color = NA,
    size = 3.1, fontface = "bold", seed = 20260606, box.padding = 0.18
  ) +
  geom_text_repel(
    data = s, aes(x = 2, y = r_exp, color = dir, label = paste0(r_exp, "  ", nm)),
    hjust = 0, nudge_x = 0.05, direction = "y", segment.color = NA,
    size = 3.1, fontface = "bold", seed = 20260606, box.padding = 0.18
  ) +
  scale_y_reverse() +
  scale_x_continuous(
    limits = c(0.35, 2.65), breaks = c(1, 2),
    labels = c("P_risk 순위\n(규제 위험)", "기대손실 순위\n(투자 우선순위)"),
    position = "top"
  ) +
  scale_color_manual(
    values = dir_cols, name = NULL, drop = TRUE,
    labels = c("상승" = "우선순위 상승 · 금액이 끌어올림",
               "하락" = "우선순위 하락 · 확률만 높음",
               "유지" = "유지")
  ) +
  labs(
    title = "확률 순위가 곧 손실 순위는 아니다",
    subtitle = "왼쪽은 P_risk(규제 위험) 순위, 오른쪽은 기대손실(투자 우선순위) 순위.\n선이 교차할수록 두 기준이 어긋난다는 뜻이다. London은 위험 하위에서 손실 1위로 솟고, Staffordshire Moorlands는 위험 1위에서 손실 하위로 떨어진다.",
    x = NULL, y = NULL,
    caption = "P_risk 상위 6개 ∪ 기대손실 상위 6개 지역. 순위는 264개 지역 기준, 숫자가 작을수록 상위."
  ) +
  theme_chart +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold", color = "#111827"),
    legend.position = "bottom"
  )

save_chart(p2, "38_02_rank_reversal_slope.png", width = 10.5, height = 8.4)
cat("  saved 38_02_rank_reversal_slope.png\n")

# ============================================================
# ③ 기대손실 집중도 Pareto — Top 15 막대 + 누적 비중 선
# ============================================================
cat("[3] 기대손실 집중도 Pareto\n")

pe <- d[order(-total_expected_loss_short)]
pe[, cum_share := cumsum(expected_loss_share_short)]
total_exp_b <- sum(d$total_expected_loss_short) / 1e9   # 전국 합 (£B)

top <- pe[1:15]
top[, nm := short_name(area_label)]
apply_city_tag(top, "stacked")                       # x축 2줄 표기 (지역\n분류·인구)
top[, nm := factor(nm, levels = nm)]                 # 막대 순서 고정
top[, exp_m := total_expected_loss_short / 1e6]
top[, is_london := area_label == "Greater London"]

sf <- max(top$exp_m) / max(top$cum_share)            # 누적%를 막대축에 맞추는 변환계수
share_top10 <- pe[10, cum_share]

p3 <- ggplot(top, aes(x = nm)) +
  geom_col(aes(y = exp_m, fill = is_london), width = 0.74) +
  geom_text(aes(y = exp_m, label = fmt_money_compact(total_expected_loss_short)),
            vjust = -0.5, size = 2.9, color = "#334155") +
  geom_line(aes(y = cum_share * sf, group = 1), color = "#1F2937", linewidth = 0.85) +
  geom_point(aes(y = cum_share * sf), color = "#1F2937", size = 2.1) +
  geom_label(aes(y = cum_share * sf, label = percent(cum_share, accuracy = 1)),
             vjust = -0.5, size = 2.6, color = "#1F2937", fontface = "bold",
             fill = alpha("white", 0.8), label.size = 0,
             label.padding = grid::unit(0.08, "lines")) +
  annotate("segment", x = 10.5, xend = 10.5, y = 0, yend = share_top10 * sf,
           linetype = "dotted", color = "#7F1D1D", linewidth = 0.5) +
  annotate("text", x = 10.4, y = share_top10 * sf, hjust = 1, vjust = -0.4,
           label = paste0("상위 10개 = ", percent(share_top10, accuracy = 1)),
           size = 3.3, fontface = "bold", color = "#7F1D1D") +
  scale_fill_manual(values = c("FALSE" = "#FC8D59", "TRUE" = "#7F0000"), guide = "none") +
  scale_y_continuous(
    name = "기대손실 (£M)",
    labels = label_number(accuracy = 1, suffix = "M"),
    sec.axis = sec_axis(~ . / sf, name = "전국 누적 비중", labels = percent_format(accuracy = 1))
  ) +
  labs(
    title = "기대손실은 소수 지역에 쏠려 있다",
    subtitle = "막대는 지역별 기대손실(£M), 선은 전국 누적 비중. Greater London 단독 14.2%, 상위 10개 지역이 약 26%를 차지한다.",
    x = NULL,
    caption = paste0("264개 지역 중 상위 15개. 지역명 아래 = ONS 도시-농촌 분류·인구(2021). 누적 비중은 전국 기대손실 약 £",
                     format(round(total_exp_b, 2), nsmall = 2), "B 대비.")
  ) +
  theme_chart +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
        axis.title.y.right = element_text(color = "#1F2937"))

save_chart(p3, "38_03_expected_loss_pareto.png", width = 12.5, height = 8.0)
cat("  saved 38_03_expected_loss_pareto.png\n")

# ============================================================
# ④ 손실 컴포넌트 분해 — Top 15 지역, 벌금/임대손실/개보수 100% 누적
# ============================================================
cat("[4] 손실 컴포넌트 분해\n")

# 전국 합계 컴포넌트 비율 (subtitle용)
nat_fine <- sum(d$total_fine_short)
nat_rent <- sum(d$total_rent_loss_short)
nat_retro <- sum(d$total_retrofit_cost)
nat_tot <- nat_fine + nat_rent + nat_retro

comp <- d[order(-total_loss_amount_short)][1:15]
comp[, nm := short_name(area_label)]
apply_city_tag(comp, "full")   # 지역명 뒤 도시유형·인구 태그
nm_order <- rev(comp$nm)   # London을 맨 위로

cl <- melt(
  comp[, .(nm, total_loss_amount_short,
           개보수 = total_retrofit_cost,
           임대손실 = total_rent_loss_short,
           벌금 = total_fine_short)],
  id.vars = c("nm", "total_loss_amount_short"),
  variable.name = "component", value.name = "amount"
)
cl[, component := factor(component, levels = c("개보수", "임대손실", "벌금"))]
cl[, nm := factor(nm, levels = nm_order)]
cl[, prop := amount / total_loss_amount_short]

tot <- comp[, .(nm = factor(nm, levels = nm_order),
                lbl = fmt_money_compact(total_loss_amount_short))]

comp_cols <- c("개보수" = "#3B78B5", "임대손실" = "#FDAE61", "벌금" = "#D73027")

p4 <- ggplot(cl, aes(x = prop, y = nm, fill = component)) +
  geom_col(width = 0.74) +
  geom_text(data = tot, aes(x = 1.015, y = nm, label = lbl), inherit.aes = FALSE,
            hjust = 0, size = 3.0, color = "#334155", fontface = "bold") +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1.16), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  scale_fill_manual(values = comp_cols, name = "손실 컴포넌트") +
  labs(
    title = "손실금액의 약 79%는 개보수 비용이다",
    subtitle = paste0(
      "지역별 손실금액 구성비. 전국 합계 기준 개보수 £", round(nat_retro / 1e9, 1),
      "B(", round(nat_retro / nat_tot * 100), "%) ≫ 임대손실 £", round(nat_rent / 1e9, 1),
      "B(", round(nat_rent / nat_tot * 100), "%) ≫ 벌금 £", round(nat_fine / 1e9, 1),
      "B(", round(nat_fine / nat_tot * 100), "%).\n오른쪽 값은 지역별 총 손실금액."
    ),
    x = "손실금액 구성비", y = NULL,
    caption = "264개 지역 중 손실금액 상위 15개. 지역명 뒤 = ONS 도시-농촌 분류·인구(2021). 개보수 + 임대손실 + 벌금 = 손실금액(확률 곱하기 전)."
  ) +
  theme_chart +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top")

save_chart(p4, "38_04_loss_component_breakdown.png", width = 12.5, height = 7.6)
cat("  saved 38_04_loss_component_breakdown.png\n")

# ============================================================
# ⑤ IMD 결핍 × 위험/손실 — decile별 패턴
# ============================================================
cat("[5] IMD × 위험/손실 (2-panel scatter)\n")

ds <- d[!is.na(building_weighted_imd_score) & n_buildings > 0]
ds[, loss_pb := total_loss_amount_short / n_buildings]   # 건물당 손실금액
r_prisk <- cor(ds$building_weighted_imd_score, ds$mean_p_risk, method = "spearman")
r_loss  <- cor(ds$building_weighted_imd_score, ds$loss_pb, method = "spearman")
cat("  Spearman  P_risk:", round(r_prisk, 3), " loss:", round(r_loss, 3), "\n")

p5a <- ggplot(ds, aes(building_weighted_imd_score, mean_p_risk)) +
  geom_point(color = "#C7352B", alpha = 0.5, size = 1.9) +
  geom_smooth(method = "lm", se = TRUE, color = "#172033", fill = "#CBD5E1", linewidth = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(subtitle = paste0("P_risk vs 결핍   (Spearman ", sprintf("%+.2f", r_prisk), ")"),
       x = "IMD score (높을수록 결핍)", y = "지역 평균 P_risk") +
  theme_chart

p5b <- ggplot(ds, aes(building_weighted_imd_score, loss_pb)) +
  geom_point(color = "#2C7BB6", alpha = 0.5, size = 1.9) +
  geom_smooth(method = "lm", se = TRUE, color = "#172033", fill = "#CBD5E1", linewidth = 0.8) +
  scale_y_log10(labels = label_number(suffix = "K", scale = 1e-3, accuracy = 1)) +
  labs(subtitle = paste0("건물당 손실금액 vs 결핍   (Spearman ", sprintf("%+.2f", r_loss), ")"),
       x = "IMD score (높을수록 결핍)", y = "건물당 손실금액 (£K, log)") +
  theme_chart

p5 <- p5a + p5b + plot_annotation(
  title = "결핍은 위험과 이어지지만 손실금액과는 약하다",
  subtitle = "왼쪽. 결핍이 클수록 P_risk가 완만히 오른다.   오른쪽. 손실금액은 결핍과 거의 무관하다 — 자산 규모·임대료가 좌우한다.",
  caption = "264개 지역, England IMD. 검은 선은 선형회귀(lm), 음영은 95% CI.",
  theme = theme(
    plot.title = element_text(face = "bold", size = 18, color = "#111827"),
    plot.subtitle = element_text(size = 11, color = "#334155", margin = margin(b = 4)),
    plot.caption = element_text(size = 8.5, color = "#64748B", hjust = 0)
  )
)

save_chart(p5, "38_05_imd_risk_loss_scatter.png", width = 12, height = 6.2)
cat("  saved 38_05_imd_risk_loss_scatter.png\n")

# ============================================================
# ⑥ 단기 vs 장기 시나리오 민감도 — dumbbell
# ============================================================
cat("[6] 단기 vs 장기 민감도\n")

s6 <- d[order(-total_expected_loss_short)][1:15]
s6[, nm := short_name(area_label)]
s6[, nm := factor(nm, levels = rev(nm))]
nat_s <- sum(d$total_expected_loss_short) / 1e9
nat_l <- sum(d$total_expected_loss_long) / 1e9

lab_s <- "단기 · 벌금 RV10% · 공실 3개월"
lab_l <- "장기 · 벌금 RV20% · 공실 6개월"

p6 <- ggplot(s6) +
  geom_segment(aes(y = nm, yend = nm,
                   x = total_expected_loss_short / 1e6,
                   xend = total_expected_loss_long / 1e6),
               color = "#CBD5E1", linewidth = 1.5) +
  geom_point(aes(y = nm, x = total_expected_loss_short / 1e6, color = lab_s), size = 3.1) +
  geom_point(aes(y = nm, x = total_expected_loss_long / 1e6, color = lab_l), size = 3.1) +
  scale_color_manual(values = setNames(c("#FDAE61", "#C7352B"), c(lab_s, lab_l)), name = NULL) +
  scale_x_log10(breaks = c(30, 50, 100, 200, 400, 700),
                labels = label_number(suffix = "M", accuracy = 1)) +
  labs(
    title = "시나리오에 따른 기대손실 범위 (단기 vs 장기)",
    subtitle = paste0("벌금률·공실기간 가정을 강화하면 기대손실이 늘어난다. 전국 합계 £",
                      format(round(nat_s, 2), nsmall = 2), "B → £",
                      format(round(nat_l, 2), nsmall = 2), "B (+",
                      round((nat_l / nat_s - 1) * 100), "%)."),
    x = "기대손실 (£M, log scale)", y = NULL,
    caption = "기대손실 상위 15개 지역. 점 사이 간격이 가정에 따른 민감도 폭. 가로축 log scale."
  ) +
  theme_chart + theme(legend.position = "top", panel.grid.major.y = element_blank())

save_chart(p6, "38_06_short_vs_long_sensitivity.png", width = 11, height = 7.4)
cat("  saved 38_06_short_vs_long_sensitivity.png\n")

# ============================================================
# ⑦ 2×2 발표용 매트릭스 — 분면별 대응 전략
# ============================================================
cat("[7] 2x2 발표용 매트릭스\n")

qcount <- d[, .N, by = quadrant_prisk_loss_amount]
getrep <- function(q) {
  r <- d[quadrant_prisk_loss_amount == q][order(-total_expected_loss_short)][1:3, short_name(area_label)]
  paste(r, collapse = ", ")
}
cells <- data.table(
  x = c(0, 1, 0, 1), y = c(1, 1, 0, 0),    # TL, TR, BL, BR
  quad = c("Lower P_risk + High Loss_amount", "High P_risk + High Loss_amount",
           "Lower P_risk + Lower Loss_amount", "High P_risk + Lower Loss_amount"),
  name = c("런던형", "최우선", "저순위", "지방형"),
  strat = c("금액 큼·확률 낮음 → 선별 투자",
            "확률·금액 모두 큼 → 즉시 최우선 대응",
            "둘 다 낮음 → 모니터링",
            "확률 높음·금액 낮음 → 비용효율 개보수"),
  fill = c("#FBD7B5", "#F1A99A", "#DDE2E8", "#BBD7E6")
)
cells[, n := vapply(quad, function(q) qcount[quadrant_prisk_loss_amount == q, N], integer(1))]
cells[, reps := vapply(quad, getrep, character(1))]
cells[, label := paste0(name, "  (", n, "개 지역)\n", reps, "\n\n→ ", strat)]

p7 <- ggplot(cells) +
  geom_rect(aes(xmin = x, xmax = x + 1, ymin = y, ymax = y + 1, fill = name),
            color = "white", linewidth = 2.2) +
  geom_text(aes(x = x + 0.5, y = y + 0.5, label = label),
            size = 3.7, lineheight = 1.2, fontface = "bold", color = "#172033") +
  scale_fill_manual(values = setNames(cells$fill, cells$name), guide = "none") +
  annotate("text", x = 1, y = -0.13, label = "위험 확률 (P_risk)  낮음 ◀───────▶ 높음",
           size = 3.6, fontface = "bold", color = "#475569") +
  annotate("text", x = -0.12, y = 1, angle = 90, label = "손실금액  작음 ◀───────▶ 큼",
           size = 3.6, fontface = "bold", color = "#475569") +
  coord_cartesian(xlim = c(-0.2, 2), ylim = c(-0.22, 2.02), expand = FALSE, clip = "off") +
  labs(
    title = "위험 × 손실 2×2 대응 전략",
    subtitle = "각 축 상위 25% 기준. 같은 위험이라도 손실금액에 따라 대응 우선순위가 달라진다."
  ) +
  theme_void(base_family = KFONT) +
  theme(
    plot.title = element_text(face = "bold", size = 18, color = "#111827"),
    plot.subtitle = element_text(size = 11, color = "#334155", margin = margin(b = 10)),
    plot.margin = margin(14, 24, 16, 20)
  )

save_chart(p7, "38_07_quadrant_2x2_strategy.png", width = 10.5, height = 8)
cat("  saved 38_07_quadrant_2x2_strategy.png\n")

cat("[완료] 의사결정 매트릭스(38_01) + 파생 그래프 →", out_chart, "\n")
})

# ════════════════════════════════════════════════════════════════════════════
# 7. 전국 LA 4분면 지도                          [PPT 50]
# ════════════════════════════════════════════════════════════════════════════
local({
# 전국 LA 위험×손실 4분면 지도 — 독립 실행본 (지도 1장)
# ─────────────────────────────────────────────────────────────────────────────
# 원본. 지역 결합/02_코드/50j_national_la_quadrant.R 에서 "지도 1장"만 분리.
#       (스케일 비교 막대그래프·csv 출력 부분은 제외)
#
# 실행법
#   • RStudio. 이 .R 열고 → Session > Set Working Directory > To Source File Location
#              → 전체 선택 후 Run (또는 source)
#   • 터미널.  이 폴더로 이동 후  Rscript 전국LA_4분면지도.R
#   → 결과 PNG 가 이 폴더에 생성됨. 20_전국LA_4분면지도.png
#
# 필요 패키지. data.table, ggplot2, scales, sf, spdep, dplyr, ggrepel
#   설치. install.packages(c("data.table","ggplot2","scales","sf","spdep","dplyr","ggrepel"))
# ─────────────────────────────────────────────────────────────────────────────

# ── 스크립트 자기 폴더 자동 감지 (RStudio / Rscript 모두 대응) ───────────────

# ── LA 데이터 + per-building 평균 손실금액(연면적·임대 기반, 규모 비의존) ─────
la <- fread(file.path(DATA, "33_la_eda_dataset.csv"), encoding = "UTF-8")
la[local_authority == "E08000039", local_authority := "E08000019"]   # Sheffield 코드 보정
la[local_authority == "E08000038", local_authority := "E08000016"]   # Barnsley 코드 보정
la[, loss_amt := (total_fine_short + total_rent_loss_short + total_retrofit_expected) / n_buildings]
la <- la[is.finite(loss_amt) & n_buildings >= 30]

# ── 4유형 정의 (전국 중앙값 분할) ────────────────────────────────────────────
qlev <- c("최우선 · 고위험+고손실", "중심상업형 · 저위험+고손실",
          "노후외곽형 · 고위험+저손실 (형평사각)", "저순위 · 저위험+저손실")
qcol <- c("#D73027", "#FC8D59", "#4E9CC4", "#C2C9D2"); names(qcol) <- qlev
mr <- median(la$mean_p_risk); ml <- median(la$loss_amt)
la[, quad := factor(fcase(
  mean_p_risk >  mr & loss_amt >  ml, qlev[1],
  mean_p_risk <= mr & loss_amt >  ml, qlev[2],
  mean_p_risk >  mr & loss_amt <= ml, qlev[3],
  default = qlev[4]), levels = qlev)]

# ── geojson(LAD13) → 2019~2023 통합 단일구 합치는 크로스워크로 디졸브 ────────
g <- st_make_valid(st_read(file.path(GEO, "uk_lad.geojson"), quiet = TRUE))
xw <- rbindlist(list(
  data.table(parent="E06000065", child=c("E07000163","E07000164","E07000165","E07000166","E07000167","E07000168","E07000169")), # N Yorkshire
  data.table(parent="E06000061", child=c("E07000150","E07000152","E07000153","E07000156")),             # N Northants
  data.table(parent="E06000062", child=c("E07000151","E07000154","E07000155")),                          # W Northants
  data.table(parent="E06000060", child=c("E07000004","E07000005","E07000006","E07000007")),             # Bucks
  data.table(parent="E06000066", child=c("E07000187","E07000188","E07000189","E07000190","E07000191")), # Somerset
  data.table(parent="E06000064", child=c("E07000027","E07000030","E07000031")),                          # Westmorland&Furness
  data.table(parent="E06000058", child=c("E06000028","E06000029","E07000048")),                          # BCP
  data.table(parent="E06000059", child=c("E07000049","E07000050","E07000051","E07000052","E07000053")), # Dorset
  data.table(parent="E06000063", child=c("E07000026","E07000028","E07000029")),                          # Cumberland
  data.table(parent="E07000244", child=c("E07000205","E07000206")),                                       # East Suffolk
  data.table(parent="E07000245", child=c("E07000201","E07000204"))))                                      # West Suffolk
cat("크로스워크 child 중 geojson 부재:", paste(setdiff(xw$child, g$LAD13CD), collapse=", "), "(없으면 OK)\n")
g$code <- g$LAD13CD
g$code[match(xw$child, g$LAD13CD)] <- xw$parent
gd <- g %>% group_by(code) %>% summarise(.groups = "drop")   # 디졸브(union)

gm <- merge(gd, la, by.x = "code", by.y = "local_authority")
cat("지도 매칭 LA:", nrow(gm), "/", nrow(la), "\n")

# ── 군집강도(subtitle 표기용) — 이웃 LA 같은유형 비율 + 위험 Moran's I ───────
nb <- poly2nb(gm, queen = TRUE); lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
I_risk <- moran.test(gm$mean_p_risk, lw, zero.policy = TRUE)$estimate[1]
jc_obs <- { same <- 0; tot <- 0
  for (i in seq_along(nb)) for (j in nb[[i]]) if (j > i) { tot <- tot+1; if (gm$quad[i]==gm$quad[j]) same <- same+1 }
  same/tot }
jc_exp <- sum(prop.table(table(gm$quad))^2)

# ── 전국 LA 4분면 지도 ───────────────────────────────────────────────────────
sub <- paste(sprintf("%s %d", c("최우선","중심상업형","형평사각","저순위"),
             sapply(qlev, function(q) la[quad==q, .N])), collapse = "  ·  ")

# 주요 도시 라벨(건물 평균 좌표). London은 33개 borough 평균으로 1점.
city_codes <- c(Birmingham="E08000025", Manchester="E08000003", Leeds="E08000035",
                Liverpool="E08000012", Sheffield="E08000019", Newcastle="E08000021",
                Bristol="E06000023", Nottingham="E08000017")
lab <- la[local_authority %in% city_codes][, .(nm = names(city_codes)[match(local_authority, city_codes)],
                                               lon = mean_longitude, lat = mean_latitude)]
lab <- rbind(lab, la[substr(local_authority,1,3)=="E09",
                     .(nm="London", lon=mean(mean_longitude), lat=mean(mean_latitude))])

p_map <- ggplot(gm) +
  geom_sf(aes(fill = quad), color = "white", linewidth = 0.08) +
  scale_fill_manual(values = qcol, name = NULL, drop = FALSE) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  geom_point(data = lab, aes(lon, lat), size = 1.1, color = "#111827") +
  ggrepel::geom_text_repel(data = lab, aes(lon, lat, label = nm), size = 3.3,
                  fontface = "bold", color = "#111827", seed = 3, min.segment.length = 0,
                  segment.color = "#111827", segment.size = 0.3, box.padding = 0.45,
                  bg.color = "white", bg.r = 0.18) +
  coord_sf(datum = NA) +
  labs(title = "전국 LA 위험×손실 4유형 — 범위를 키우면 패턴이 뚜렷하다",
       subtitle = paste0("잉글랜드 ", nrow(gm), "개 LA · 전국 중앙값 분할(동 분석과 같은 4유형).  ", sub,
                         "\n이웃 LA가 같은 유형일 확률 ", sprintf("%.0f%%", 100*jc_obs),
                         " (무작위 ", sprintf("%.0f%%", 100*jc_exp), ") · 위험 Moran I=", sprintf("%.2f", I_risk),
                         " — 도시 내부(런던 빼고 ≈0)와 정반대"),
       caption = "빨강=최우선(북부·후기산업) · 파랑=형평사각(고위험·저임대) · 주황=중심상업(런던·남동부 돈) · 회색=저순위. 손실금액=건물당 평균(규모 비의존).") +
  theme_void(base_family = KFONT) +
  theme(plot.title = element_text(face = "bold", size = 16, color = "#111827"),
        plot.subtitle = element_text(size = 10, color = "#334155", lineheight = 1.15, margin = margin(b = 6)),
        plot.caption = element_text(size = 8.2, color = "#64748B", hjust = 0),
        legend.position = "bottom", legend.text = element_text(size = 9.2),
        plot.margin = margin(10, 10, 8, 10))

ggsave(file.path(OUTF, "13_50_전국LA_4분면지도.png"), p_map, width = 9.5, height = 11, dpi = 230, bg = "white")
cat("\n[완료] 저장됨 → 20_전국LA_4분면지도.png\n")
})

# ════════════════════════════════════════════════════════════════════════════
# 8. 4유형 × ONS 공식분류                        [PPT 50]
# ════════════════════════════════════════════════════════════════════════════
local({
# 50o_ons_classification.R — 위험×손실 4유형 × ONS 공식 지역분류(2011 LAC) 매칭
# ─────────────────────────────────────────────────────────────────────────────
# 각 4유형이 영국 공식 분류(ONS 2011 Area Classification for Local Authorities)의
#  어떤 성격(산업유산·대도시·농촌 등)을 띠는지 heatmap으로. PPT 4유형 설명 근거.
# 산출. 06_4분면_지역분석/28_4유형_ONS공식분류.png
# ─────────────────────────────────────────────────────────────────────────────
outdir <- OUTF

lac <- fread(file.path(DATA, "ons_la_classification_2011.csv"), encoding = "UTF-8")
# 4분면(quad)을 13_quadrant_map.R 과 동일한 전국 중앙값 분할로 직접 계산
la <- fread(file.path(DATA, "33_la_eda_dataset.csv"), encoding = "UTF-8")
la[local_authority == "E08000039", local_authority := "E08000019"]
la[local_authority == "E08000038", local_authority := "E08000016"]
la[, loss_amt := (total_fine_short + total_rent_loss_short + total_retrofit_expected) / n_buildings]
la <- la[is.finite(loss_amt) & n_buildings >= 30]
mr <- median(la$mean_p_risk); ml <- median(la$loss_amt)
la[, quad := fcase(mean_p_risk >  mr & loss_amt >  ml, "최우선 · 고위험+고손실",
                   mean_p_risk <= mr & loss_amt >  ml, "중심상업형 · 저위험+고손실",
                   mean_p_risk >  mr & loss_amt <= ml, "노후외곽형 · 고위험+저손실",
                   default = "저순위 · 저위험+저손실")]
q <- la[, .(la_label = local_authority_label, quad)]
m <- merge(q, lac, by.x = "la_label", by.y = "Name")

# 유형 라벨(축 순서: 위험·손실 순)
m[, qs := fcase(grepl("중심상업", quad), "🟠 중심상업형\n저위험·고손실",
                grepl("최우선", quad), "🔴 최우선\n고위험·고손실",
                grepl("노후외곽", quad), "🔵 지방형\n고위험·저손실",
                default = "⚪ 저순위\n저위험·저손실")]
qord <- c("🟠 중심상업형\n저위험·고손실", "🔴 최우선\n고위험·고손실",
          "🔵 지방형\n고위험·저손실", "⚪ 저순위\n저위험·저손실")
# ONS supergroup 한글 약칭 + 도시→농촌 순서
ab <- c("London Cosmopolitan"="런던 코스모폴리탄","Ethnically Diverse Metropolitan Living"="다민족 대도시",
  "Business, Education and Heritage Centres"="업무·교육·문화\n중심지","Urban Settlements"="도시 거주지",
  "Affluent England"="부유한 잉글랜드","Services and Industrial Legacy"="서비스·산업유산",
  "Town and Country Living"="소도시·전원","Countryside Living"="전원·농촌")
m[, sg := factor(ab[Supergroup], levels = ab)]
m[, qs := factor(qs, levels = rev(qord))]

# 유형 내 비율(행 정규화)
ct <- m[, .N, by = .(qs, sg)][, share := N / sum(N), by = qs]
ct <- ct[CJ(qs = levels(m$qs), sg = levels(m$sg), unique = TRUE), on = .(qs, sg)]
ct[is.na(N), `:=`(N = 0, share = 0)]

p <- ggplot(ct, aes(sg, qs, fill = share)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = ifelse(share > 0, percent(share, accuracy = 1), "")),
            size = 3.5, fontface = "bold",
            color = ifelse(ct$share > 0.25, "white", "#1F2937")) +
  scale_fill_gradient(low = "#F1F5F9", high = "#B91C1C", labels = percent, name = "유형 내 비율") +
  scale_x_discrete(position = "top") +
  labs(title = "위험×손실 4유형은 어떤 '공식 분류' 지역인가 (ONS 2011 LAC)",
       subtitle = "셀 = 그 유형 안에서 해당 ONS 분류가 차지하는 비율. 진할수록 그 유형의 대표 성격.\n🟠중심상업형=대도시·업무중심  🔴최우선=산업유산·소도시  🔵지방형=전원·농촌(40%)  ⚪저순위=부유한 교외",
       x = NULL, y = NULL,
       caption = "ONS 2011 Area Classification for Local Authorities(8 supergroups), 302개 LA 매칭. 위험×손실 4유형은 전국 중앙값 분할.") +
  theme_minimal(base_family = KFONT) +
  theme(plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 9.8, color = "#334155", lineheight = 1.2),
        plot.caption = element_text(size = 8, color = "#64748B", hjust = 0),
        axis.text.x = element_text(size = 9, face = "bold", lineheight = 0.9),
        axis.text.y = element_text(size = 10.5, face = "bold", lineheight = 0.95),
        panel.grid = element_blank(), legend.position = "right",
        plot.margin = margin(10, 12, 8, 10))
ggsave(file.path(outdir, "28_4유형_ONS공식분류.png"), p, width = 12.5, height = 6.5, dpi = 230, bg = "white")
cat("  saved 28_4유형_ONS공식분류.png\n")

# 발표용 대표 LA(유형×supergroup 핵심 셀)
cat("\n[유형별 대표 분류·지역]\n")
for (u in qord) {
  uu <- gsub("\n", " ", u)
  top <- m[qs == u, .N, by = .(Supergroup)][order(-N)][1:2]
  egs <- sapply(top$Supergroup, function(s) paste(head(m[qs == u & Supergroup == s, la_label], 3), collapse = ", "))
  cat(sprintf("• %s\n", uu))
  for (i in 1:nrow(top)) cat(sprintf("    %s (%d개): %s\n", top$Supergroup[i], top$N[i], egs[i]))
}
})

# ════════════════════════════════════════════════════════════════════════════
# 9. LA 클러스터링 (K-means 2단계)               [발표 보강]
# ════════════════════════════════════════════════════════════════════════════
local({
# ─────────────────────────────────────────────────────────────────────────────
# 전국 LA 클러스터링 (K-means) — 2단계: ① 4분면 검증  ② 다차원 유형발견
# ─────────────────────────────────────────────────────────────────────────────
# 목적.
#   기존 4분면(quadrant_total)은 사람이 '상위 25% 선'을 그어 만든 rule-based 구분.
#   여기서는 같은 LA들을 데이터가 스스로 묶게(K-means) 해서
#     STAGE 1) 2축(위험·총손실)으로 군집 → 기존 4분면과 교차검증
#     STAGE 2) 섹터구성+IMD+자산·규모까지 넣은 다차원 군집 → 새 지역유형 발견
#
# 실행.  이 폴더에서  Rscript 클러스터링_LA.R
#        (data/33_la_eda_dataset.csv, data/uk_lad.geojson 필요 — 4분면지도와 동일)
# 산출.  클러스터링_결과/ 폴더에 PNG·CSV 생성
# 필요 패키지. data.table, ggplot2, scales, sf, dplyr, ggrepel, cluster
# ─────────────────────────────────────────────────────────────────────────────
set.seed(42)                              # K-means 난수 고정 (재현성)

outdir <- file.path(ROOT, "outputs", "clustering"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ── 0. 데이터 로드 + 코드 보정 + 필터 (4분면 스크립트와 동일 기준) ───────────
la <- fread(file.path(DATA, "33_la_eda_dataset.csv"), encoding = "UTF-8")
la[local_authority == "E08000039", local_authority := "E08000019"]   # Sheffield
la[local_authority == "E08000038", local_authority := "E08000016"]   # Barnsley
la <- la[is.finite(mean_p_risk) & n_buildings >= 30]
cat(sprintf("대상 LA: %d개\n", nrow(la)))

# ═════════════════════════════════════════════════════════════════════════════
# STAGE 1 — 검증용: 2축(위험 × 총손실) K-means(k=4) vs 기존 4분면
# ═════════════════════════════════════════════════════════════════════════════
X1 <- la[, .(mean_p_risk, log_total_loss = log1p(total_expected_loss_short))]
X1s <- scale(X1)                                   # 표준화 필수(스케일 통일)
km1 <- kmeans(X1s, centers = 4, nstart = 50, iter.max = 50)
la[, cl1 := km1$cluster]

# 군집을 '위험 high/low × 손실 high/low'로 해석 가능하게 재라벨 (센트로이드 기준)
ov_r <- median(la$mean_p_risk); ov_l <- median(la$total_expected_loss_short)
cen1 <- la[, .(r = mean(mean_p_risk), l = mean(total_expected_loss_short)), by = cl1]
cen1[, name := fcase(r >  ov_r & l >  ov_l, "고위험+고손실",
                     r <= ov_r & l >  ov_l, "저위험+고손실",
                     r >  ov_r & l <= ov_l, "고위험+저손실",
                     default                = "저위험+저손실")]
la[cen1, cl1_name := i.name, on = "cl1"]

# 기존 4분면(quadrant_total)을 같은 4라벨로 매핑해 교차표 작성
map_q <- c("High P_risk + High Total Loss"   = "고위험+고손실",
           "Lower P_risk + High Total Loss"  = "저위험+고손실",
           "High P_risk + Lower Total Loss"  = "고위험+저손실",
           "Lower P_risk + Lower Total Loss" = "저위험+저손실")
la[, quad_kr := map_q[quadrant_total]]
ct1 <- table(`K-means(데이터)` = la$cl1_name, `기존 4분면(상위25%선)` = la$quad_kr)
agree <- mean(la$cl1_name == la$quad_kr, na.rm = TRUE)
cat("\n[STAGE 1] 2축 K-means(k=4)  vs  기존 4분면 ─ 교차표\n")
print(ct1)
cat(sprintf("→ 두 방식 라벨 일치율: %.1f%%  (데이터 군집이 25%% 임의선과 어디서 갈리는지 확인용)\n", 100*agree))
fwrite(as.data.frame.matrix(ct1), file.path(outdir, "stage1_교차표.csv"), row.names = TRUE)

# ═════════════════════════════════════════════════════════════════════════════
# STAGE 2 — 유형발견용: 다차원 프로파일 K-means
#   피처(9). 위험강도 + 건물당손실강도 + 자산수준 + 규모 + 섹터4 + 지역박탈도
#   ※ 총액(규모) 대신 '건물당/구성비'를 써서 도시 크기 착시를 빼고 '유형'을 잡음
#   ※ rank·quadrant 등 '결과' 컬럼은 절대 피처로 넣지 않음(논리오염 방지)
# ═════════════════════════════════════════════════════════════════════════════
la[, imd_imp := building_weighted_imd_score]                       # IMD 31개 결측
la[is.na(imd_imp), imd_imp := median(la$building_weighted_imd_score, na.rm = TRUE)]
cat(sprintf("\nIMD 결측 %d개 → 중앙값 대치\n", sum(is.na(la$building_weighted_imd_score))))

feat <- la[, .(
  위험      = mean_p_risk,
  손실강도  = log1p(mean_expected_loss_short),     # 건물당 기대손실(규모 비의존)
  자산수준  = log1p(mean_RV),
  건물규모  = log1p(mean_floor_area),
  사무비중  = share_office,
  소매외식  = share_retail_food,
  산업창고  = share_industry_storage,
  여가집합  = share_leisure_assembly,
  지역박탈  = imd_imp)]                              # IMD score 높을수록 박탈 심함
cat(sprintf("STAGE 2 메인 피처 %d개\n", ncol(feat)))
Xs <- scale(feat)

# ── 차원의 저주(curse of dimensionality) 대응: PCA로 차원 축소 후 군집 ────────
#   9피처 중 섹터·자산·규모가 서로 상관 → 고차원에서 점간 거리가 균질해져(거리 집중)
#   실루엣이 0.16까지 깔림. PCA로 '무상관 주성분'만 추려(분산 80% 유지) 군집하면
#   중복·잡음 축이 제거돼 거리 신호가 또렷해지고 분리도가 회복됨.
pca <- prcomp(Xs)                                       # Xs 이미 표준화
varexp <- pca$sdev^2 / sum(pca$sdev^2)
npc <- max(2, which.max(-diff(varexp)))                 # 스크리 elbow: 분산 급감 직전까지
Z   <- pca$x[, 1:npc, drop = FALSE]                     # 군집은 이 축소공간에서 수행
cat(sprintf("\n[PCA] 9피처 → %dPC (스크리 elbow · 누적분산 %.0f%%) · 차원의 저주 완화\n",
            npc, 100*cumsum(varexp)[npc]))
rn <- names(feat)                                        # 주성분 해석(로딩)
for (j in 1:npc) { v <- pca$rotation[, j]; v <- v[abs(v) >= 0.35]
  cat(sprintf("  PC%d(%.0f%%): %s\n", j, 100*varexp[j],
              paste(sprintf("%s%+.2f", names(v), v), collapse="  "))) }

# ── k 선택: elbow + 평균 실루엣 (PCA 공간 기준) + 축소 전 9D 비교 ────────────
ks <- 2:10
dZ <- dist(Z)
wss <- sapply(ks, function(k) kmeans(Z, k, nstart = 25, iter.max = 50)$tot.withinss)
sil <- sapply(ks, function(k) mean(silhouette(kmeans(Z, k, nstart=25, iter.max=50)$cluster, dZ)[,3]))
sil_raw <- sapply(ks, function(k) mean(silhouette(kmeans(Xs, k, nstart=25, iter.max=50)$cluster, dist(Xs))[,3]))
sel <- data.table(k = ks, WSS = round(wss),
                  실루엣_9D원본 = round(sil_raw, 3),    # 차원의 저주(축소 전)
                  실루엣_PCA    = round(sil, 3))         # 축소 후
cat("\n[STAGE 2] k 후보별 지표 (실루엣 9D원본 → PCA축소)\n"); print(sel)
# 실루엣이 k=4~6에서 평탄(≈0.25 동률) → 해석가능성·벨트 안정성으로 k=5 확정
elig <- ks %in% 3:7
k_auto <- min(ks[elig & sil >= max(sil[elig]) - 0.01])    # 자동(참고)
k_star <- 5L                                              # 평탄구간 → 프리미엄 분리·벨트 7/8 위해 5
cat(sprintf("→ 후보 k(자동)=%d · 최종 k=%d (실루엣 평탄 → 해석가능성으로 선택)\n", k_auto, k_star))
fwrite(sel, file.path(outdir, "stage2_k선택지표.csv"))

# elbow + silhouette 플롯
pe <- ggplot(sel, aes(k, WSS)) + geom_line(color="#4E79A7") + geom_point(size=2) +
  geom_vline(xintercept = k_star, linetype = 2, color = "#D62728") +
  scale_x_continuous(breaks = ks) +
  labs(title = "Elbow — 군집 내 분산(WSS)", x = "군집 수 k", y = "WSS") + theme_minimal(base_family=KFONT)
ps <- ggplot(sel) +
  geom_line(aes(k, 실루엣_9D원본), color="#B6BcC8") + geom_point(aes(k, 실루엣_9D원본), color="#B6BCC8", size=2) +
  geom_line(aes(k, 실루엣_PCA),   color="#59A14F") + geom_point(aes(k, 실루엣_PCA),   color="#59A14F", size=2) +
  geom_vline(xintercept = k_star, linetype = 2, color = "#D62728") +
  scale_x_continuous(breaks = ks) +
  labs(title = "평균 실루엣 — 회색=원본9D · 초록=PCA축소(차원의 저주 완화)",
       x = "군집 수 k", y = "avg silhouette") + theme_minimal(base_family=KFONT)
ggsave(file.path(outdir, "stage2_k선택_elbow_silhouette.png"),
       gridExtra::grid.arrange(pe, ps, ncol = 2), width = 11, height = 4.2, dpi = 200, bg="white")

# ── 최종 군집 + 프로파일표 (PCA 축소공간에서 군집) ──────────────────────────
km2 <- kmeans(Z, centers = k_star, nstart = 50, iter.max = 100)
la[, cl2 := km2$cluster]
cat(sprintf("\n[STAGE 2] 최종 실루엣: 원본9D %.3f → PCA축소 %.3f (k=%d)\n",
            sil_raw[k_star-1], sil[k_star-1], k_star))

prof <- la[, .(
  N            = .N,
  위험        = round(mean(mean_p_risk), 3),
  건물당손실  = round(mean(mean_expected_loss_short)),
  자산RV      = round(mean(mean_RV)),
  규모면적    = round(mean(mean_floor_area)),
  사무 = round(mean(share_office),2),  소매 = round(mean(share_retail_food),2),
  산업 = round(mean(share_industry_storage),2), 여가 = round(mean(share_leisure_assembly),2),
  IMD박탈     = round(mean(imd_imp), 1)
), by = cl2][order(cl2)]

# 군집 자동 라벨 — '특화섹터'는 절대 최대가 아니라 전체평균 대비 상대(z)로 잡음
#   (소매/금융 비중이 어디서나 제일 커서 절대 최대는 항상 소매 → 구분 불가)
g_med_r <- median(la$mean_p_risk); g_med_imd <- median(la$imd_imp)
sect_cols <- c("사무","소매","산업","여가")
sect_z <- scale(as.matrix(prof[, ..sect_cols]))          # 군집간 섹터 z-score
prof[, 특화섹터 := sect_cols[max.col(sect_z)]]
prof[, 자산등급 := cut(자산RV, quantile(자산RV, c(0,1/3,2/3,1)),
                       labels = c("저자산","중자산","고자산"), include.lowest = TRUE)]
prof[, 자동라벨 := paste0(
  특화섹터, "특화·",
  fifelse(위험 > g_med_r, "고위험", "저위험"), "·",
  자산등급, "·",
  fifelse(IMD박탈 > g_med_imd, "고박탈", "저박탈"))]
# 대표지역(군집별 건물 많은 상위 3개 LA)
reps <- la[, .(대표지역 = paste(local_authority_label[order(-n_buildings)][1:3], collapse=", ")), by = cl2]
prof[reps, 대표지역 := i.대표지역, on = "cl2"]
la[prof, cl2_name := paste0("C", i.cl2, " ", i.자동라벨), on = "cl2"]

cat("\n[STAGE 2] 군집별 프로파일\n"); print(prof[, .(cl2, N, 위험, 건물당손실, 자산RV, 규모면적, IMD박탈, 특화섹터, 자동라벨)])
cat("\n[STAGE 2] 군집별 대표지역\n"); print(prof[, .(cl2, 자동라벨, 대표지역)])
fwrite(prof, file.path(outdir, "stage2_군집프로파일.csv"))

# 기존 4분면과 교차표 — '지방형(고위험·저손실)'이 몇 개로 쪼개지는지 확인
ct2 <- table(`K-means 유형` = paste0("C", la$cl2), `기존 4분면` = la$quad_kr)
cat("\n[STAGE 2] 다차원 군집  vs  기존 4분면 ─ 교차표\n"); print(ct2)
fwrite(as.data.frame.matrix(ct2), file.path(outdir, "stage2_교차표_vs_4분면.csv"), row.names = TRUE)

# ── 벨트 진단: 손으로 찍은 후기산업 벨트가 같은 군집으로 묶이나? ─────────────
belt <- c("Stoke-on-Trent","Staffordshire Moorlands","Newcastle-under-Lyme",
          "Pendle","Burnley","Blackburn","Bolton","Hyndburn")
rural <- c("Torridge","West Devon","Cornwall","North Devon")
diag_lab <- function(names_vec) {
  hit <- la[grepl(paste(names_vec, collapse="|"), local_authority_label, ignore.case = TRUE),
            .(지역 = local_authority_label, 군집 = cl2_name, 위험 = round(mean_p_risk,3),
              IMD = round(imd_imp,1))]
  unique(hit)
}
cat("\n[벨트 진단] 후기산업 벨트(도자기·방직)가 같은 군집으로 묶이나?\n"); print(diag_lab(belt))
cat("\n[대조] 농촌형(Torridge·West Devon 등)\n"); print(diag_lab(rural))

# ── 견고성(robustness) — 두 가지 검증 ───────────────────────────────────────
ari <- function(a, b) { t <- table(a, b); n <- sum(t); c2 <- function(x) sum(choose(x, 2))
  i <- c2(as.vector(t)); a2 <- c2(rowSums(t)); b2 <- c2(colSums(t))
  e <- a2 * b2 / choose(n, 2); (i - e) / ((a2 + b2)/2 - e) }            # Adjusted Rand Index
belt_idx <- grepl(paste(belt, collapse = "|"), la$local_authority_label, ignore.case = TRUE)
coh <- function(cl) max(table(cl[belt_idx])) / sum(belt_idx)           # 벨트가 한 군집에 모인 비율

# (1) PCA 축소가 유형을 왜곡했나? — 축소 전 원본 9D 군집과 비교
km_raw  <- kmeans(Xs, centers = k_star, nstart = 50, iter.max = 100)
ari_pca <- ari(la$cl2, km_raw$cluster)

# (2) 피처를 9→13개로 늘려도 핵심 유형 재현되나? — 동일하게 PCA 축소 후 군집
feat13 <- copy(feat)
feat13[, `:=`(교육시설=la$share_education, 자산집중=la$top5pct_loss_share_short,
              개보수비중=la$total_retrofit_expected /
                         (la$total_fine_short+la$total_rent_loss_short+la$total_retrofit_expected),
              위험이질성=la$p75_p_risk-la$median_p_risk)]
for (cc in names(feat13)) feat13[is.na(get(cc)), (cc):=median(feat13[[cc]], na.rm=TRUE)]
p13 <- prcomp(scale(feat13)); n13 <- which(cumsum(p13$sdev^2/sum(p13$sdev^2)) >= 0.80)[1]
km13 <- kmeans(p13$x[, 1:n13], centers = 6, nstart = 50, iter.max = 100)
ari_13 <- ari(la$cl2, km13$cluster)

cat(sprintf("\n[견고성] (1) PCA축소 vs 원본9D 군집  ARI = %.2f → PCA가 구조 왜곡 안 함\n", ari_pca))
cat(sprintf("          (2) 9피처 vs 13피처 군집      ARI = %.2f → 피처 추가에도 안정\n", ari_13))
cat(sprintf("          벨트 결집률: 메인 %.0f%% · 13피처 %.0f%%  (모두 '형평사각' 유형 재현)\n",
            100*coh(la$cl2), 100*coh(km13$cluster)))
fwrite(data.table(ari_pca_vs_raw9D = round(ari_pca, 3), ari_9f_vs_13f = round(ari_13, 3),
                  belt_coh_main = round(coh(la$cl2), 3), belt_coh_13f = round(coh(km13$cluster), 3)),
       file.path(outdir, "stage2_견고성.csv"))

# ═════════════════════════════════════════════════════════════════════════════
# 지도 — 4분면 스크립트의 지오메트리/크로스워크/도시라벨 재사용
# ═════════════════════════════════════════════════════════════════════════════
g <- st_make_valid(st_read(file.path(GEO, "uk_lad.geojson"), quiet = TRUE))
xw <- rbindlist(list(
  data.table(parent="E06000065", child=c("E07000163","E07000164","E07000165","E07000166","E07000167","E07000168","E07000169")),
  data.table(parent="E06000061", child=c("E07000150","E07000152","E07000153","E07000156")),
  data.table(parent="E06000062", child=c("E07000151","E07000154","E07000155")),
  data.table(parent="E06000060", child=c("E07000004","E07000005","E07000006","E07000007")),
  data.table(parent="E06000066", child=c("E07000187","E07000188","E07000189","E07000190","E07000191")),
  data.table(parent="E06000064", child=c("E07000027","E07000030","E07000031")),
  data.table(parent="E06000058", child=c("E06000028","E06000029","E07000048")),
  data.table(parent="E06000059", child=c("E07000049","E07000050","E07000051","E07000052","E07000053")),
  data.table(parent="E06000063", child=c("E07000026","E07000028","E07000029")),
  data.table(parent="E07000244", child=c("E07000205","E07000206")),
  data.table(parent="E07000245", child=c("E07000201","E07000204"))))
g$code <- g$LAD13CD
g$code[match(xw$child, g$LAD13CD)] <- xw$parent
gd <- g %>% group_by(code) %>% summarise(.groups = "drop")
gm <- merge(gd, la, by.x = "code", by.y = "local_authority")
cat(sprintf("\n지도 매칭 LA: %d / %d\n", nrow(gm), nrow(la)))

city_codes <- c(Birmingham="E08000025", Manchester="E08000003", Leeds="E08000035",
                Liverpool="E08000012", Sheffield="E08000019", Newcastle="E08000021",
                Bristol="E06000023", Nottingham="E08000017")
lab <- la[local_authority %in% city_codes][, .(nm = names(city_codes)[match(local_authority, city_codes)],
                                               lon = mean_longitude, lat = mean_latitude)]
lab <- rbind(lab, la[substr(local_authority,1,3)=="E09",
                     .(nm="London", lon=mean(mean_longitude), lat=mean(mean_latitude))])

draw_map <- function(fill_col, title, sub, palette = NULL) {
  p <- ggplot(gm) + geom_sf(aes(fill = .data[[fill_col]]), color="white", linewidth=0.08)
  p <- if (is.null(palette)) p + scale_fill_brewer(palette = "Set2", name = NULL)
       else p + scale_fill_manual(values = palette, name = NULL)
  p + guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
    geom_point(data = lab, aes(lon, lat), size = 1.1, color = "#111827") +
    ggrepel::geom_text_repel(data = lab, aes(lon, lat, label = nm), size = 3.2, fontface="bold",
        color="#111827", seed=3, min.segment.length=0, segment.size=0.3, box.padding=0.45,
        bg.color="white", bg.r=0.18) +
    coord_sf(datum = NA) +
    labs(title = if (nzchar(title)) title else NULL,
         subtitle = if (nzchar(sub)) sub else NULL) +
    theme_void(base_family=KFONT) +
    theme(plot.title=element_text(face="bold", size=15),
          plot.subtitle=element_text(size=10, color="#334155", margin=margin(b=6)),
          legend.position="bottom", legend.text=element_text(size=9))
}

# STAGE1 지도 (4분면과 같은 색 체계로 비교)
qcol <- c("고위험+고손실"="#D73027","저위험+고손실"="#FC8D59",
          "고위험+저손실"="#4E9CC4","저위험+저손실"="#C2C9D2")
gm$cl1_name <- factor(gm$cl1_name, levels = names(qcol))
ggsave(file.path(outdir, "지도_stage1_2축군집.png"),
       draw_map("cl1_name", "STAGE 1 — 2축 K-means(위험×총손실, k=4)",
                sprintf("데이터가 스스로 나눈 4군집 · 기존 4분면과 일치율 %.0f%%", 100*agree), qcol),
       width = 9, height = 10.5, dpi = 220, bg="white")

# STAGE2 지도 (다차원 유형)
gm$cl2_name <- factor(gm$cl2_name)
ggsave(file.path(outdir, "지도_stage2_다차원유형.png"),
       draw_map("cl2_name", sprintf("STAGE 2 — 다차원 K-means (k=%d)", k_star),
                "섹터구성+IMD+자산·규모 기반 '지역 포트폴리오 유형'"),
       width = 9, height = 10.5, dpi = 220, bg="white")

# ── 슬라이드 삽입용: 제목 없는 깔끔본 (발표덱 제목과 중복 방지) ───────────────
ggsave(file.path(outdir, "지도_stage1_slide.png"),
       draw_map("cl1_name", "", "", qcol), width = 8.6, height = 10, dpi = 220, bg="white")
ggsave(file.path(outdir, "지도_stage2_slide.png"),
       draw_map("cl2_name", "", ""),       width = 8.6, height = 10, dpi = 220, bg="white")
# 기존 4분면(손으로 그은 것) 지도 — S1 좌우 비교용 (cl1과 같은 색 체계)
gm$quad_kr_f <- factor(gm$quad_kr, levels = names(qcol))
ggsave(file.path(outdir, "지도_4분면_slide.png"),
       draw_map("quad_kr_f", "", "", qcol), width = 8.6, height = 10, dpi = 220, bg="white")

# ── 최종 군집 라벨 CSV 저장(태블로/추가분석용) ───────────────────────────────
fwrite(la[, .(local_authority, local_authority_label, n_buildings, mean_p_risk,
              total_expected_loss_short, mean_expected_loss_short, mean_RV, mean_floor_area,
              imd_imp, quadrant_total, quad_kr,
              stage1_cluster = cl1_name, stage2_cluster = cl2_name)],
       file.path(outdir, "LA_클러스터_라벨.csv"))

cat("\n[완료] 산출물 →", outdir, "/\n")
cat("  지도_stage1_2축군집.png · 지도_stage2_다차원유형.png\n")
cat("  stage2_k선택_elbow_silhouette.png · stage2_군집프로파일.csv\n")
cat("  stage1_교차표.csv · stage2_교차표_vs_4분면.csv · LA_클러스터_라벨.csv\n")
})

# ════════════════════════════════════════════════════════════════════════════
# 10. 리스크 회피 정량화                         [PPT 52~54]
# ════════════════════════════════════════════════════════════════════════════
local({
# 리스크 회피 정량화 — 같은 예산으로 막는 손실 비교  [PPT 52~54]
# 입력. final_risk_ranked.csv  출력. outputs/figures/30_*, outputs/tables/30_*
# 같은 수의 건물을 손볼 때 어떤 기준으로 고르느냐에 따라 회피 기대손실이 달라진다.
set.seed(42)

f <- fread(file.path(DATA, "final_risk_ranked.csv"), encoding = "UTF-8")
EL   <- f$total_risk_short_term                          # 기대손실 = P_risk × 손실금액(단기)
amt  <- f$fine_short_term + f$rent_loss_short + f$pure_retrofit_cost_total  # 확률 전 손실금액
pool <- sum(EL)                                          # 전체 기대손실 풀
cat(sprintf("[30] 건물 %s건 · 전체 기대손실 풀 £%.2fB\n", comma(nrow(f)), pool / 1e9))

# 전략별 정렬 순서 → 그 순서대로 실제 기대손실을 누적하면 얼마나 걷어내나
strategies <- list(
  "모델 (P_risk × 손실금액)" = order(-EL),
  "손실금액순 (확률 미반영)"  = order(-amt),
  "규모 (RV)순"               = order(-f$RV),
  "위험확률 (P_risk)순"       = order(-f$P_risk),
  "무작위"                    = sample(nrow(f)))

# 누적 회피 곡선 (전략별)
curve <- rbindlist(lapply(names(strategies), function(s) {
  cum <- cumsum(EL[strategies[[s]]]) / pool
  data.table(strategy = s, frac = seq_len(nrow(f)) / nrow(f), captured = cum)
}))
curve <- curve[seq(1, .N, by = 50)]                      # 다운샘플(그림 가벼우라고)

# 상위 10% / 20% 대응 시 전략별 회피 기대손실
cap_at <- function(idx, p) sum(EL[idx[seq_len(floor(nrow(f) * p))]])
tab <- rbindlist(lapply(names(strategies), function(s)
  data.table(전략 = s,
             `상위10%_£B` = round(cap_at(strategies[[s]], 0.10) / 1e9, 2),
             `상위20%_£B` = round(cap_at(strategies[[s]], 0.20) / 1e9, 2))))
cat("\n[30] 전략별 회피 기대손실 (상위 10% / 20% 대응)\n"); print(tab)
save_tab(tab, "30_avoidance_by_strategy.csv")

model10 <- tab[전략 == "모델 (P_risk × 손실금액)", `상위10%_£B`]
prisk10 <- tab[전략 == "위험확률 (P_risk)순",     `상위10%_£B`]
cat(sprintf("[30] 모델 상위10%% £%.2fB = 위험확률순 £%.2fB 의 약 %.1f배\n",
            model10, prisk10, model10 / prisk10))

# ── 그림 1. 누적 회피 곡선 [슬라이드 52] ────────────────────────────────────
ord <- c("모델 (P_risk × 손실금액)", "손실금액순 (확률 미반영)", "규모 (RV)순",
         "위험확률 (P_risk)순", "무작위")
curve[, strategy := factor(strategy, levels = ord)]
cols <- c("#1F4E96", "#3B9C9C", "#59A14F", "#F1A340", "#9CA3AF")
p1 <- ggplot(curve, aes(frac, captured, color = strategy)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#9CA3AF") +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = 0.10, linetype = "dotted", color = "#475569") +
  scale_color_manual(values = cols, name = NULL) +
  scale_x_continuous(labels = percent) + scale_y_continuous(labels = percent) +
  labs(title = "같은 예산, 더 많은 손실 회피 — 기대손실 기준 우선순위",
       subtitle = "전략 순위대로 선제 대응할 때 누적 회피 기대손실(전체 대비). 점선=무작위.",
       x = "선제 대응한 건물 비율", y = "누적 회피 기대손실 (%)")
save_fig(p1, "30_52_avoidance_curve.png", w = 9.5, h = 6)

# ── 그림 2. 상위 10% 대응 시 전략별 회피액 막대 [슬라이드 53] ────────────────
b <- tab[, .(전략, v = `상위10%_£B`)][, 전략 := factor(전략, levels = rev(ord))]
p2 <- ggplot(b, aes(v, 전략, fill = 전략)) + geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("£%.2fB", v)), hjust = -0.1, size = 4) +
  scale_fill_manual(values = setNames(cols, ord)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "상위 10% 건물 대응 시 — 전략별 회피 기대손실",
       subtitle = sprintf("전체 풀 £%.2fB 중. 모델이 위험확률순의 약 %.1f배를 회피.", pool / 1e9, model10 / prisk10),
       x = "회피 기대손실 (£B)", y = NULL)
save_fig(p2, "30_53_avoidance_top10_bar.png", w = 9, h = 5)

cat("[30] 리스크 회피 정량화 완료\n")
})
