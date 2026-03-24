pkgs <- c(
  "haven", "dplyr", "tidyr", "purrr", "ggplot2", "splines",
  "nnet", "sandwich", "lmtest", "boot", "broom", "scales",
  "readr", "flextable", "officer", "patchwork", "logistf",
  "marginaleffects"
)

to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}

invisible(lapply(pkgs, library, character.only = TRUE))

options(scipen = 999)
set.seed(42)
theme_set(theme_minimal(base_size = 11))

out_dir <- "outputs_koweps_final"
if(!dir.exists(out_dir)) dir.create(out_dir)

# ============================================================
# [1] 데이터 로드
# ============================================================

h15 <- read_dta('/Users/eunseongpark/Downloads/한국복지패널_1_19차_차수별_데이터(release250327)stata.zip/(2020년 15차 한국복지패널조사) 데이터(beta5)_stata/koweps_h15_2020_beta5.dta')
h16 <- read_dta('/Users/eunseongpark/Downloads/한국복지패널_1_19차_차수별_데이터(release250327)stata.zip/(2021년 16차 한국복지패널조사) 데이터(beta4)_stata/koweps_h16_2021_beta4.dta')
h17 <- read_dta('/Users/eunseongpark/Downloads/한국복지패널_1_19차_차수별_데이터(release250327)stata.zip/(2022년 17차 한국복지패널조사) 데이터(beta3)_stata/koweps_h17_2022_beta3.dta')
h18 <- read_dta('/Users/eunseongpark/Downloads/한국복지패널_1_19차_차수별_데이터(release250327)stata.zip/(2023년 18차 한국복지패널조사) 데이터(beta2)_stata/koweps_h18_2023_beta2.dta')
h19 <- read_dta('/Users/eunseongpark/Downloads/한국복지패널_1_19차_차수별_데이터(release250327)stata.zip/(2024년 19차 한국복지패널조사) 데이터(beta1)_stata/koweps_h19_2024_beta1.dta')

cat("✓ 5개 연도 데이터 로드 완료\n\n")

# ============================================================
# [2] 기본 변수 생성 함수
# ============================================================
make_base <- function(df, wave, year){
  id_col   <- paste0("h", wave, "_merkey")
  inc_col  <- paste0("h", wave, "_din")
  prem_col <- paste0("h", wave, "05_3aq2")
  ins_cols <- paste0("h", wave, "02_3aq", 1:9)
  
  needed <- c(id_col, inc_col, prem_col, ins_cols)
  miss <- setdiff(needed, names(df))
  if(length(miss) > 0) stop("Wave ", wave, " missing: ", paste(miss, collapse=", "))
  
  df %>%
    transmute(
      hhid = .data[[id_col]],
      ins  = rowSums(across(all_of(ins_cols)), na.rm = TRUE),
      inc  = .data[[inc_col]],
      prm  = .data[[prem_col]] * 12,
      brd  = prm / inc,
      year = year
    )
}

# ============================================================
# [3] 완전패널 식별
# ============================================================
ids_all <- Reduce(
  intersect,
  list(
    h15$h15_merkey, h16$h16_merkey, h17$h17_merkey,
    h18$h18_merkey, h19$h19_merkey
  )
)

cat(sprintf("✓ 완전패널 가구: %d개\n\n", length(ids_all)))

# ============================================================
# [4] 기본 데이터 생성
# ============================================================
b15 <- make_base(h15, 15, 2020) %>% filter(hhid %in% ids_all)
b16 <- make_base(h16, 16, 2021) %>% filter(hhid %in% ids_all)
b17 <- make_base(h17, 17, 2022) %>% filter(hhid %in% ids_all)
b18 <- make_base(h18, 18, 2023) %>% filter(hhid %in% ids_all)
b19 <- make_base(h19, 19, 2024) %>% filter(hhid %in% ids_all)

# ============================================================
# [5] 전환 데이터 생성
# ============================================================
make_transition <- function(bt, bt1){
  inner_join(bt, bt1, by = "hhid", suffix = c("_t", "_t1")) %>%
    transmute(
      hhid     = hhid,
      year_t   = year_t,
      ins_t    = ins_t,
      inc_t    = inc_t,
      prm_t    = prm_t,
      brd_t    = brd_t,
      ins_t1   = ins_t1,
      decrease = as.integer(ins_t1 < ins_t),
      exit     = as.integer(ins_t >= 1 & ins_t1 == 0)
    )
}

tr_20_21 <- make_transition(b15, b16)
tr_21_22 <- make_transition(b16, b17)
tr_22_23 <- make_transition(b17, b18)
tr_23_24 <- make_transition(b18, b19)

transitions_all <- bind_rows(tr_20_21, tr_21_22, tr_22_23, tr_23_24)

cat(sprintf("✓ 총 전환 관측치: %d개\n", nrow(transitions_all)))
cat(sprintf("  - 2020→2021: %d\n", nrow(tr_20_21)))
cat(sprintf("  - 2021→2022: %d\n", nrow(tr_21_22)))
cat(sprintf("  - 2022→2023: %d\n", nrow(tr_22_23)))
cat(sprintf("  - 2023→2024: %d\n\n", nrow(tr_23_24)))

# ============================================================
# [6] 소득그룹 분류 (20/60/20)
# ============================================================
transitions_all <- transitions_all %>%
  group_by(year_t) %>%
  mutate(
    inc_rank  = rank(inc_t, ties.method = "random"),
    n_t       = n(),
    inc_group = case_when(
      inc_rank <= 0.2 * n_t  ~ "Low20",
      inc_rank <= 0.8 * n_t  ~ "Mid60",
      TRUE                   ~ "High20"
    )
  ) %>%
  ungroup() %>%
  select(-inc_rank, -n_t)

# ============================================================
# [7] 분석표본 고정 (필터 적용 금지)
# ============================================================
analysis_sample_all <- transitions_all %>%
  filter(
    ins_t >= 1,
    inc_t > 0,
    is.finite(brd_t)
  ) %>%
  mutate(
    prm_zero    = as.integer(prm_t == 0),
    year_f      = factor(year_t),
    inc_group_f = factor(inc_group, levels = c("Low20", "Mid60", "High20")),
    action      = case_when(
      exit == 1 ~ "Exit",
      decrease == 1 & exit == 0 ~ "DecreaseOnly",
      TRUE ~ "NoChange"
    ),
    action = factor(action, levels = c("NoChange", "DecreaseOnly", "Exit"))
  )

cat(sprintf("✓ 분석표본 확정: %d개 관측치\n", nrow(analysis_sample_all)))
cat(sprintf("  - NoChange: %d (%.1f%%)\n",
            sum(analysis_sample_all$action == "NoChange"),
            100 * mean(analysis_sample_all$action == "NoChange")))
cat(sprintf("  - DecreaseOnly: %d (%.1f%%)\n",
            sum(analysis_sample_all$action == "DecreaseOnly"),
            100 * mean(analysis_sample_all$action == "DecreaseOnly")))
cat(sprintf("  - Exit: %d (%.1f%%)\n",
            sum(analysis_sample_all$action == "Exit"),
            100 * mean(analysis_sample_all$action == "Exit")))
cat(sprintf("  - prm_zero share: %.1f%%\n\n",
            100 * mean(analysis_sample_all$prm_zero)))

# DecreaseOnly 부분표본
sub_decOnly <- analysis_sample_all %>%
  filter(exit == 0) %>%
  mutate(dec_only = as.integer(decrease == 1))

# ============================================================
# [8] 기술통계 계산
# ============================================================
tech_overall <- analysis_sample_all %>%
  summarise(
    N               = n(),
    decrease_rate   = mean(decrease),
    exit_rate       = mean(exit),
    mean_brd        = mean(brd_t),
    p95_brd         = quantile(brd_t, 0.95),
    p99_brd         = quantile(brd_t, 0.99),
    prm_zero_share  = mean(prm_zero)
  )

tech_by_group <- analysis_sample_all %>%
  group_by(inc_group) %>%
  summarise(
    N = n(),
    decrease_rate = mean(decrease),
    exit_rate = mean(exit),
    prm_zero_share = mean(prm_zero),
    .groups = "drop"
  )

tech_by_year <- analysis_sample_all %>%
  group_by(year_t) %>%
  summarise(
    N = n(),
    decrease_rate = mean(decrease),
    exit_rate = mean(exit),
    mean_brd = mean(brd_t),
    p95_brd = quantile(brd_t, 0.95),
    p99_brd = quantile(brd_t, 0.99),
    prm_zero_share = mean(prm_zero),
    .groups = "drop"
  )

# CSV 저장
write_csv(tech_overall, file.path(out_dir, "T_tech_overall.csv"))
write_csv(tech_by_group, file.path(out_dir, "T_tech_by_group.csv"))
write_csv(tech_by_year, file.path(out_dir, "T_tech_by_year.csv"))

cat("✓ 기술통계 계산 완료\n\n")

# ============================================================
# [9] 메인 모형: 다항로짓 (MNL)
# ============================================================
m_mnl <- nnet::multinom(
  action ~ brd_t + ins_t + inc_group_f + year_f + prm_zero,
  data = analysis_sample_all,
  trace = FALSE
)

mnl_sum <- summary(m_mnl)
mnl_coef <- as.matrix(mnl_sum$coefficients)
mnl_se <- as.matrix(mnl_sum$standard.errors)
storage.mode(mnl_coef) <- "double"
storage.mode(mnl_se) <- "double"

mnl_z <- mnl_coef / mnl_se
mnl_p <- 2 * (1 - pnorm(abs(mnl_z)))

mnl_tbl <- data.frame()
for(i in 1:nrow(mnl_coef)) {
  for(j in 1:ncol(mnl_coef)) {
    mnl_tbl <- rbind(mnl_tbl, data.frame(
      outcome = colnames(mnl_coef)[j],
      term = rownames(mnl_coef)[i],
      coef = mnl_coef[i, j],
      se = mnl_se[i, j],
      z = mnl_z[i, j],
      p = mnl_p[i, j],
      OR = exp(mnl_coef[i, j])
    ))
  }
}

write_csv(mnl_tbl, file.path(out_dir, "T_MNL_results.csv"))
cat("✓ MNL 모형 적합 완료\n\n")

# ============================================================
# [10] 보조 모형: 이항로짓 (클러스터 강건)
# ============================================================

# Exit 모형
m_exit_bin <- glm(
  exit ~ brd_t + ins_t + inc_group_f + year_f + prm_zero,
  family = binomial,
  data = analysis_sample_all
)

V_exit <- vcovCL(m_exit_bin, cluster = analysis_sample_all$hhid)
exit_ct <- coeftest(m_exit_bin, vcov. = V_exit)
exit_tbl <- data.frame(
  term = rownames(exit_ct),
  estimate = exit_ct[, 1],
  se = exit_ct[, 2],
  z = exit_ct[, 3],
  p = exit_ct[, 4],
  OR = exp(exit_ct[, 1])
)

write_csv(exit_tbl, file.path(out_dir, "T_Exit_logit.csv"))

# DecreaseOnly 모형
m_dec_bin <- glm(
  dec_only ~ brd_t + ins_t + inc_group_f + year_f + prm_zero,
  family = binomial,
  data = sub_decOnly
)

V_dec <- vcovCL(m_dec_bin, cluster = sub_decOnly$hhid)
dec_ct <- coeftest(m_dec_bin, vcov. = V_dec)
dec_tbl <- data.frame(
  term = rownames(dec_ct),
  estimate = dec_ct[, 1],
  se = dec_ct[, 2],
  z = dec_ct[, 3],
  p = dec_ct[, 4],
  OR = exp(dec_ct[, 1])
)

write_csv(dec_tbl, file.path(out_dir, "T_DecreaseOnly_logit.csv"))

cat("✓ 이항로짓 모형 적합 완료\n\n")

# ============================================================
# [11] 비선형성 검정 (스플라인)
# ============================================================

# Exit 선형 vs 스플라인
m_exit_spl <- glm(
  exit ~ ns(brd_t, 3) + ins_t + inc_group_f + year_f + prm_zero,
  family = binomial,
  data = analysis_sample_all
)

lrt_exit <- anova(m_exit_bin, m_exit_spl, test = "Chisq")

# DecreaseOnly 선형 vs 스플라인
m_dec_spl <- glm(
  dec_only ~ ns(brd_t, 3) + ins_t + inc_group_f + year_f + prm_zero,
  family = binomial,
  data = sub_decOnly
)

lrt_dec <- anova(m_dec_bin, m_dec_spl, test = "Chisq")

cat(sprintf("✓ Exit 비선형성 검정 (LRT p=%.4f)\n", lrt_exit$`Pr(>Chi)`[2]))
cat(sprintf("✓ DecreaseOnly 비선형성 검정 (LRT p=%.4f)\n\n", lrt_dec$`Pr(>Chi)`[2]))

# ============================================================
# [12] 임계점 자동탐색 (2단계)
# ============================================================

tau_stage1 <- seq(0.0005, 0.05, by = 0.0005)
tau_stage2 <- seq(0.05, 0.30, by = 0.005)
tau_grid <- sort(unique(c(tau_stage1, tau_stage2)))

scan_tau_exit_precise <- function(data) {
  groups <- unique(data$inc_group)
  out <- list()
  
  for(g in groups) {
    sub <- data %>% filter(inc_group == g)
    
    res_g <- map_dfr(tau_grid, function(tau) {
      sub$split <- as.integer(sub$brd_t > tau)
      
      fit <- try(
        glm(exit ~ brd_t + split + ins_t + factor(year_t),
            family = binomial, data = sub),
        silent = TRUE
      )
      
      if(inherits(fit, "try-error") || isFALSE(fit$converged)) {
        return(data.frame(group = g, tau = tau, aic = Inf))
      }
      data.frame(group = g, tau = tau, aic = AIC(fit))
    })
    
    out[[g]] <- res_g
  }
  
  bind_rows(out)
}

scan_tau_dec_precise <- function(data) {
  groups <- unique(data$inc_group)
  out <- list()
  
  for(g in groups) {
    sub <- data %>% filter(inc_group == g)
    
    res_g <- map_dfr(tau_grid, function(tau) {
      sub$split <- as.integer(sub$brd_t > tau)
      
      fit <- try(
        glm(dec_only ~ brd_t + split + ins_t + factor(year_t),
            family = binomial, data = sub),
        silent = TRUE
      )
      
      if(inherits(fit, "try-error") || isFALSE(fit$converged)) {
        return(data.frame(group = g, tau = tau, aic = Inf))
      }
      data.frame(group = g, tau = tau, aic = AIC(fit))
    })
    
    out[[g]] <- res_g
  }
  
  bind_rows(out)
}

# Exit 임계점 탐색
aic_exit2 <- scan_tau_exit_precise(analysis_sample_all)
best_tau_exit2 <- aic_exit2 %>%
  group_by(group) %>%
  slice_min(aic, n = 1, with_ties = FALSE) %>%
  ungroup()

# DecreaseOnly 임계점 탐색
aic_dec2 <- scan_tau_dec_precise(sub_decOnly)
best_tau_dec2 <- aic_dec2 %>%
  group_by(group) %>%
  slice_min(aic, n = 1, with_ties = FALSE) %>%
  ungroup()

write_csv(best_tau_exit2, file.path(out_dir, "T_best_tau_exit.csv"))
write_csv(best_tau_dec2, file.path(out_dir, "T_best_tau_decreaseonly.csv"))

cat("✓ Exit 임계점 탐색 완료:\n")
print(best_tau_exit2)
cat("\n✓ DecreaseOnly 임계점 탐색 완료:\n")
print(best_tau_dec2)
cat("\n")

# ============================================================
# FIGURE 1: 전환 행동 분포 (전체)
# ============================================================

fig1_data <- analysis_sample_all %>%
  mutate(outcome3 = case_when(
    exit == 1 ~ "Exit",
    exit == 0 & decrease == 1 ~ "DecreaseOnly",
    TRUE ~ "NoChange"
  )) %>%
  count(outcome3) %>%
  mutate(
    pct = 100 * n / sum(n),
    outcome3 = factor(outcome3,
                      levels = c("NoChange", "DecreaseOnly", "Exit"))
  )

fig1 <- ggplot(fig1_data, aes(x = outcome3, y = pct, fill = outcome3)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct, n)),
            vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(
    values = c("NoChange" = "#4A90D9", "DecreaseOnly" = "#F5A623", "Exit" = "#D0021B")
  ) +
  scale_y_continuous(limits = c(0, max(fig1_data$pct) * 1.2)) +
  labs(
    title = "Figure 1. Overall Transition Distribution",
    x = "",
    y = "Percentage (%)"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 10)
  )

ggsave(file.path(out_dir, "Fig1_overall_distribution.png"),
       fig1, width = 8, height = 6, dpi = 300)
print(fig1)

cat("✓ Figure 1 저장\n\n")

# ============================================================
# FIGURE 2: 전환 행동 분포 (소득그룹별)
# ============================================================

fig2_data <- analysis_sample_all %>%
  mutate(outcome3 = case_when(
    exit == 1 ~ "Exit",
    exit == 0 & decrease == 1 ~ "DecreaseOnly",
    TRUE ~ "NoChange"
  )) %>%
  count(inc_group, outcome3) %>%
  group_by(inc_group) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    inc_group = factor(inc_group, levels = c("Low20", "Mid60", "High20")),
    outcome3 = factor(outcome3,
                      levels = c("NoChange", "DecreaseOnly", "Exit"))
  )

fig2 <- ggplot(fig2_data, aes(x = outcome3, y = pct, fill = outcome3)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  facet_wrap(~inc_group, nrow = 1) +
  scale_fill_manual(
    values = c("NoChange" = "#4A90D9", "DecreaseOnly" = "#F5A623", "Exit" = "#D0021B")
  ) +
  scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) +
  labs(
    title = "Figure 2. Transition Distribution by Income Group",
    x = "",
    y = "Percentage (%)"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave(file.path(out_dir, "Fig2_by_income_group.png"),
       fig2, width = 12, height = 5, dpi = 300)
print(fig2)

cat("✓ Figure 2 저장\n\n")

# ============================================================
# FIGURE 3: Exit 사건수 (소득그룹×연도)
# ============================================================

fig3_data <- analysis_sample_all %>%
  filter(exit == 1) %>%
  count(inc_group, year_t) %>%
  mutate(inc_group = factor(inc_group, levels = c("Low20", "Mid60", "High20")))

fig3 <- ggplot(fig3_data, aes(x = factor(year_t), y = n, fill = inc_group)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(
    name = "Income Group",
    values = c("Low20" = "#D0021B", "Mid60" = "#F5A623", "High20" = "#4A90D9")
  ) +
  labs(
    title = "Figure 3. Exit Events by Income Group and Year",
    x = "Year",
    y = "Number of Exit Events"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "Fig3_exit_by_year.png"),
       fig3, width = 10, height = 6, dpi = 300)
print(fig3)

cat("✓ Figure 3 저장\n\n")

# ============================================================
# FIGURE 4: 보험료 부담률 분포 (소득그룹별)
# ============================================================

fig4_data <- analysis_sample_all %>%
  filter(brd_t >= 0, brd_t <= 0.30) %>%
  mutate(inc_group = factor(inc_group, levels = c("Low20", "Mid60", "High20")))

fig4 <- ggplot(fig4_data, aes(x = brd_t, fill = inc_group)) +
  geom_density(alpha = 0.65, color = NA) +
  facet_wrap(~inc_group, ncol = 1) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.30)) +
  scale_fill_manual(
    name = "Income Group",
    values = c("Low20" = "#D0021B", "Mid60" = "#F5A623", "High20" = "#4A90D9")
  ) +
  labs(
    title = "Figure 4. Premium Burden Rate Distribution by Income Group (0-30%)",
    x = "Premium Burden Rate",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(out_dir, "Fig4_burden_distribution.png"),
       fig4, width = 8, height = 8, dpi = 300)
print(fig4)

cat("✓ Figure 4 저장\n\n")

# ============================================================
# FIGURE 5: Exit 오즈비 (클러스터 강건)
# ============================================================

or_exit_df <- exit_tbl %>%
  filter(term != "(Intercept)") %>%
  mutate(
    lower = exp(estimate - 1.96 * se),
    upper = exp(estimate + 1.96 * se),
    term = reorder(term, OR)
  )

fig5 <- ggplot(or_exit_df, aes(x = OR, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray60", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.3, color = "#D0021B", linewidth = 0.9) +
  geom_point(size = 4, color = "#D0021B") +
  scale_x_log10() +
  labs(
    title = "Figure 5. Exit: Odds Ratios with 95% CI (Cluster-Robust)",
    x = "Odds Ratio (log scale)",
    y = ""
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold")
  )

ggsave(file.path(out_dir, "Fig5_exit_OR.png"),
       fig5, width = 9, height = 6, dpi = 300)
print(fig5)

cat("✓ Figure 5 저장\n\n")

# ============================================================
# FIGURE 6: DecreaseOnly 오즈비 (클러스터 강건)
# ============================================================

or_dec_df <- dec_tbl %>%
  filter(term != "(Intercept)") %>%
  mutate(
    lower = exp(estimate - 1.96 * se),
    upper = exp(estimate + 1.96 * se),
    term = reorder(term, OR)
  )

fig6 <- ggplot(or_dec_df, aes(x = OR, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray60", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = lower, xmax = upper),
                 height = 0.3, color = "#F5A623", linewidth = 0.9) +
  geom_point(size = 4, color = "#F5A623") +
  scale_x_log10() +
  labs(
    title = "Figure 6. DecreaseOnly: Odds Ratios with 95% CI (Cluster-Robust)",
    x = "Odds Ratio (log scale)",
    y = ""
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold")
  )

ggsave(file.path(out_dir, "Fig6_decreaseonly_OR.png"),
       fig6, width = 9, height = 6, dpi = 300)
print(fig6)

cat("✓ Figure 6 저장\n\n")

# ============================================================
# FIGURE 7: Exit 임계점 AIC 곡선 (전체 범위)
# ============================================================

fig7 <- ggplot(aic_exit2, aes(x = tau * 100, y = aic)) +
  geom_line(linewidth = 1) +
  facet_wrap(~group, scales = "free_y", nrow = 1) +
  geom_vline(data = best_tau_exit2, aes(xintercept = tau * 100),
             linetype = "dashed", color = "#D0021B", linewidth = 1) +
  geom_point(data = best_tau_exit2, aes(x = tau * 100, y = aic),
             size = 3, color = "#D0021B") +
  labs(
    title = "Figure 7. Exit: AIC Scan for Threshold (Full Range)",
    x = "Threshold τ (%)",
    y = "AIC"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave(file.path(out_dir, "Fig7_exit_AIC_full.png"),
       fig7, width = 12, height = 5, dpi = 300)
print(fig7)

cat("✓ Figure 7 저장\n\n")

# ============================================================
# FIGURE 8: Exit 임계점 AIC 곡선 (0-5% 확대)
# ============================================================

fig8 <- ggplot(aic_exit2 %>% filter(tau <= 0.05),
               aes(x = tau * 100, y = aic)) +
  geom_line(linewidth = 1) +
  facet_wrap(~group, scales = "free_y", nrow = 1) +
  geom_vline(data = best_tau_exit2, aes(xintercept = tau * 100),
             linetype = "dashed", color = "#D0021B", linewidth = 1) +
  geom_point(data = best_tau_exit2 %>% filter(tau <= 0.05),
             aes(x = tau * 100, y = aic),
             size = 3, color = "#D0021B") +
  labs(
    title = "Figure 8. Exit: AIC Scan for Threshold (0-5% Zoom)",
    x = "Threshold τ (%)",
    y = "AIC"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave(file.path(out_dir, "Fig8_exit_AIC_zoom.png"),
       fig8, width = 12, height = 5, dpi = 300)
print(fig8)

cat("✓ Figure 8 저장\n\n")

# ============================================================
# FIGURE 9: DecreaseOnly 임계점 AIC 곡선 (전체 범위)
# ============================================================

fig9 <- ggplot(aic_dec2, aes(x = tau * 100, y = aic)) +
  geom_line(linewidth = 1) +
  facet_wrap(~group, scales = "free_y", nrow = 1) +
  geom_vline(data = best_tau_dec2, aes(xintercept = tau * 100),
             linetype = "dashed", color = "#F5A623", linewidth = 1) +
  geom_point(data = best_tau_dec2, aes(x = tau * 100, y = aic),
             size = 3, color = "#F5A623") +
  labs(
    title = "Figure 9. DecreaseOnly: AIC Scan for Threshold (Full Range)",
    x = "Threshold τ (%)",
    y = "AIC"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave(file.path(out_dir, "Fig9_dec_AIC_full.png"),
       fig9, width = 12, height = 5, dpi = 300)
print(fig9)

cat("✓ Figure 9 저장\n\n")

# ============================================================
# FIGURE 10: DecreaseOnly 임계점 AIC 곡선 (0-5% 확대)
# ============================================================

fig10 <- ggplot(aic_dec2 %>% filter(tau <= 0.05),
                aes(x = tau * 100, y = aic)) +
  geom_line(linewidth = 1) +
  facet_wrap(~group, scales = "free_y", nrow = 1) +
  geom_vline(data = best_tau_dec2, aes(xintercept = tau * 100),
             linetype = "dashed", color = "#F5A623", linewidth = 1) +
  geom_point(data = best_tau_dec2 %>% filter(tau <= 0.05),
             aes(x = tau * 100, y = aic),
             size = 3, color = "#F5A623") +
  labs(
    title = "Figure 10. DecreaseOnly: AIC Scan for Threshold (0-5% Zoom)",
    x = "Threshold τ (%)",
    y = "AIC"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    strip.text = element_text(size = 10, face = "bold")
  )

ggsave(file.path(out_dir, "Fig10_dec_AIC_zoom.png"),
       fig10, width = 12, height = 5, dpi = 300)
print(fig10)

cat("✓ Figure 10 저장\n\n")

# ============================================================
# FIGURE 11: Exit 스플라인 예측확률 (소득그룹별)
# ============================================================

# 예측 그리드
pred_grid <- expand.grid(
  brd_t = seq(0, 0.30, length.out = 250),
  inc_group = c("Low20", "Mid60", "High20"),
  year_t = 2020
)

ins_medians <- analysis_sample_all %>%
  group_by(inc_group) %>%
  summarise(ins_med = median(ins_t), .groups = "drop")

pred_grid <- pred_grid %>%
  left_join(ins_medians, by = "inc_group") %>%
  mutate(
    ins_t = ins_med,
    prm_zero = 0,
    year_f = factor(2020, levels = levels(analysis_sample_all$year_f)),
    inc_group_f = factor(inc_group, levels = c("Low20", "Mid60", "High20"))
  )

pred_grid$pred_exit <- predict(m_exit_spl,
                               newdata = pred_grid,
                               type = "response")

fig11 <- ggplot(pred_grid, aes(x = brd_t, y = pred_exit, color = inc_group)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(limits = c(0, 0.30),
                     breaks = seq(0, 0.30, by = 0.05),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 0.35),
                     breaks = seq(0, 0.35, by = 0.05),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(
    name = "Income Group",
    values = c("Low20" = "#D0021B", "Mid60" = "#F5A623", "High20" = "#4A90D9")
  ) +
  labs(
    title = "Figure 11. Exit: Predicted Probability by Income Group (Spline)",
    x = "Premium Burden Rate",
    y = "Predicted Probability of Exit"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "Fig11_exit_spline.png"),
       fig11, width = 10, height = 6, dpi = 300)
print(fig11)

cat("✓ Figure 11 저장\n\n")

# ============================================================
# FIGURE 12: DecreaseOnly 스플라인 예측확률 (소득그룹별)
# ============================================================

pred_grid$pred_dec <- predict(m_dec_spl,
                              newdata = pred_grid %>%
                                select(-ins_med, -pred_exit),
                              type = "response")

fig12 <- ggplot(pred_grid, aes(x = brd_t, y = pred_dec, color = inc_group)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(limits = c(0, 0.30),
                     breaks = seq(0, 0.30, by = 0.05),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 0.35),
                     breaks = seq(0, 0.35, by = 0.05),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_color_manual(
    name = "Income Group",
    values = c("Low20" = "#D0021B", "Mid60" = "#F5A623", "High20" = "#4A90D9")
  ) +
  labs(
    title = "Figure 12. DecreaseOnly: Predicted Probability by Income Group (Spline)",
    x = "Premium Burden Rate",
    y = "Predicted Probability of DecreaseOnly"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

ggsave(file.path(out_dir, "Fig12_decreaseonly_spline.png"),
       fig12, width = 10, height = 6, dpi = 300)
print(fig12)

cat("✓ Figure 12 저장\n\n")

# ============================================================
# FIGURE 13: Low20 Exit 확대 (저부담 구간)
# ============================================================

fig13 <- pred_grid %>%
  filter(inc_group == "Low20", brd_t <= 0.15) %>%
  ggplot(aes(x = brd_t, y = pred_exit)) +
  geom_line(color = "#D0021B", linewidth = 1.4) +
  scale_x_continuous(limits = c(0, 0.15),
                     breaks = seq(0, 0.15, by = 0.025),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Figure 13. Exit: Low20 Low-Burden Region (0-15%)",
    x = "Premium Burden Rate",
    y = "Predicted Probability of Exit"
  ) +
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11, face = "bold")
  )

ggsave(file.path(out_dir, "Fig13_exit_low20_zoom.png"),
       fig13, width = 9, height = 6, dpi = 300)
print(fig13)

cat("✓ Figure 13 저장\n\n")

# ============================================================
# 임계점 검증: 전후 사건수/사건률
# ============================================================

exit_counts <- analysis_sample_all %>%
  left_join(best_tau_exit2, by = c("inc_group" = "group")) %>%
  mutate(above = brd_t > tau) %>%
  group_by(inc_group, above) %>%
  summarise(
    N = n(),
    exit_events = sum(exit == 1),
    exit_rate = mean(exit == 1),
    .groups = "drop"
  )

dec_counts <- sub_decOnly %>%
  left_join(best_tau_dec2, by = c("inc_group" = "group")) %>%
  mutate(above = brd_t > tau) %>%
  group_by(inc_group, above) %>%
  summarise(
    N = n(),
    dec_events = sum(dec_only == 1),
    dec_rate = mean(dec_only == 1),
    .groups = "drop"
  )

write_csv(exit_counts, file.path(out_dir, "T_exit_counts_around_threshold.csv"))
write_csv(dec_counts, file.path(out_dir, "T_decreaseonly_counts_around_threshold.csv"))

cat("✓ 임계점 검증 데이터 저장\n\n")

# ============================================================
# 최종 요약 리포트
# ============================================================

report_path <- file.path(out_dir, "FINAL_REPORT.txt")
sink(report_path, split = TRUE)

cat("========================================\n")
cat("KoWePS 보험 해지행동 분석 - 최종 보고\n")
cat("========================================\n\n")

cat(sprintf("[표본 정보]\n"))
cat(sprintf("- 완전패널 가구: %d개\n", length(ids_all)))
cat(sprintf("- 전체 전환 관측치: %d개\n", nrow(transitions_all)))
cat(sprintf("- 분석표본: %d개\n", nrow(analysis_sample_all)))
cat(sprintf("- NoChange: %d (%.1f%%)\n",
            sum(analysis_sample_all$action == "NoChange"),
            100 * mean(analysis_sample_all$action == "NoChange")))
cat(sprintf("- DecreaseOnly: %d (%.1f%%)\n",
            sum(analysis_sample_all$action == "DecreaseOnly"),
            100 * mean(analysis_sample_all$action == "DecreaseOnly")))
cat(sprintf("- Exit: %d (%.1f%%)\n",
            sum(analysis_sample_all$action == "Exit"),
            100 * mean(analysis_sample_all$action == "Exit")))
cat(sprintf("- 보험료 0인 비중: %.1f%%\n\n",
            100 * mean(analysis_sample_all$prm_zero)))

cat("[기술통계]\n")
cat("전체:\n")
print(tech_overall)
cat("\n소득그룹별:\n")
print(tech_by_group)
cat("\n연도별:\n")
print(tech_by_year)

cat("\n[비선형성 검정 (LRT)]\n")
cat(sprintf("- Exit: χ²(df=2) = %.4f, p = %.4f\n",
            lrt_exit$Deviance[2], lrt_exit$`Pr(>Chi)`[2]))
cat(sprintf("- DecreaseOnly: χ²(df=2) = %.4f, p = %.4f\n\n",
            lrt_dec$Deviance[2], lrt_dec$`Pr(>Chi)`[2]))

cat("[최적 임계점 (AIC 기준)]\n")
cat("Exit:\n")
print(best_tau_exit2)
cat("\nDecreaseOnly:\n")
print(best_tau_dec2)

cat("\n[임계점 검증: Exit]\n")
print(exit_counts)

cat("\n[임계점 검증: DecreaseOnly]\n")
print(dec_counts)

cat("\n========================================\n")
cat("분석 완료\n")
cat("========================================\n")

sink()

cat("✓ 최종 리포트 저장 완료\n\n")

# ============================================================
# 마무리
# ============================================================

cat("========================================\n")
cat("모든 분석이 완료되었습니다!\n")
cat("========================================\n\n")

cat(sprintf("출력 디렉토리: %s\n\n", out_dir))

cat("생성된 파일 목록:\n")
cat("  [그래프]\n")
cat("  - Fig1_overall_distribution.png\n")
cat("  - Fig2_by_income_group.png\n")
cat("  - Fig3_exit_by_year.png\n")
cat("  - Fig4_burden_distribution.png\n")
cat("  - Fig5_exit_OR.png\n")
cat("  - Fig6_decreaseonly_OR.png\n")
cat("  - Fig7_exit_AIC_full.png\n")
cat("  - Fig8_exit_AIC_zoom.png\n")
cat("  - Fig9_dec_AIC_full.png\n")
cat("  - Fig10_dec_AIC_zoom.png\n")
cat("  - Fig11_exit_spline.png\n")
cat("  - Fig12_decreaseonly_spline.png\n")
cat("  - Fig13_exit_low20_zoom.png\n\n")

cat("  [표]\n")
cat("  - T_tech_overall.csv\n")
cat("  - T_tech_by_group.csv\n")
cat("  - T_tech_by_year.csv\n")
cat("  - T_MNL_results.csv\n")
cat("  - T_Exit_logit.csv\n")
cat("  - T_DecreaseOnly_logit.csv\n")
cat("  - T_best_tau_exit.csv\n")
cat("  - T_best_tau_decreaseonly.csv\n")
cat("  - T_exit_counts_around_threshold.csv\n")
cat("  - T_decreaseonly_counts_around_threshold.csv\n")
cat("  - FINAL_REPORT.txt\n\n")
