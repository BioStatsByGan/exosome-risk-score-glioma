## 发表级绘图配置（舰长图规：全英文 / 无图内标题 / 刻度向内 / 四边细线 / Okabe-Ito / 矢量 PDF 可编辑）
suppressMessages({
  library(ggplot2); library(grid); library(patchwork)
})
OI <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7F7F7F", "#000000")
OI2 <- c(High = "#D55E00", Low = "#0072B2", Up = "#D55E00", Down = "#0072B2", NS = "grey78")

theme_pub <- function(base_size = 11, base_family = "sans", border = TRUE) {
  th <- theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
      panel.border = if (border) element_rect(fill = NA, colour = "black", linewidth = 0.45) else element_blank(),
      axis.line = if (border) element_blank() else element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(colour = "black", linewidth = 0.45),
      axis.ticks.length = unit(-0.12, "cm"),                 # 刻度向内
      axis.text.x = element_text(margin = margin(t = 6), colour = "black"),
      axis.text.y = element_text(margin = margin(r = 6), colour = "black"),
      axis.text = element_text(colour = "black"),
      plot.title = element_blank(),                          # 图内无标题
      plot.subtitle = element_blank(),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.title = element_text(size = base_size - 0.5),
      strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.35),
      strip.text = element_text(colour = "black"),
      plot.margin = margin(6, 8, 6, 6)
    )
  th
}
theme_set(theme_pub())

save_pdf <- function(plot, file, width, height) {
  grDevices::cairo_pdf(file, width = width, height = height, family = "sans")   # 文字保持可编辑矢量
  print(plot); dev.off()
  png_file <- sub("\\.pdf$", ".png", file)
  grDevices::png(png_file, width = width, height = height, units = "in", res = 200, type = "cairo")
  print(plot); dev.off()
  cat("saved:", file, "\n")
}

## 火山图：渐变配色（按 FDR 深度着色），标注 Top 基因
volcano_pub <- function(d, fc = "logFC", fdr = "fdr", gene = "gene", fc_cut = 1, fdr_cut = 0.05,
                        label_n = 12, xlab = expression(log[2]~"fold change"), ylab = expression(-log[10]~"FDR")) {
  d <- as.data.frame(d)
  d$.fc <- as.numeric(d[[fc]]); d$.fdr <- as.numeric(d[[fdr]])
  d$sig <- ifelse(d$.fdr < fdr_cut & d$.fc > fc_cut, "Up",
                  ifelse(d$.fdr < fdr_cut & d$.fc < -fc_cut, "Down", "NS"))
  d$score <- -log10(d$.fdr)
  sig <- d[d$sig != "NS", ]
  lab <- head(sig[order(-abs(sig$.fc)), ], label_n)
  ggplot(d, aes(.fc, score)) +
    geom_point(data = d[d$sig == "NS", ], colour = "grey85", size = 0.7, alpha = 0.8) +
    geom_point(data = d[d$sig != "NS", ], aes(colour = score, size = score), alpha = 0.95) +
    scale_colour_gradientn(colours = c("#2c7bb6", "#abd9e9", "#fdae61", "#d7191c"), name = expression(-log[10]~"FDR")) +
    scale_size_continuous(range = c(0.9, 2.8), guide = "none") +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", linewidth = 0.25, colour = "grey40") +
    geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", linewidth = 0.25, colour = "grey40") +
    labs(x = xlab, y = ylab) +
    ggrepel::geom_text_repel(data = lab, aes(label = .data[[gene]]), size = 2.6, max.overlaps = 30,
                             segment.size = 0.25, segment.colour = "grey55", min.segment.length = 0)
}

## 非对称 U 形 GSEA 曲线（替代 GO/KEGG 气泡图）
gsea_curve <- function(ranked_genes, gene_set, title = NULL, n_perm = 0, col = "#0072B2") {
  ok <- intersect(gene_set, ranked_genes)
  hit <- which(ranked_genes %in% ok)
  N <- length(ranked_genes); k <- length(ok)
  es <- cumsum(ifelse(seq_len(N) %in% hit, 1 / k, -1 / (N - k)))
  df <- data.frame(rank = seq_len(N), es = es)
  ymax <- max(es); ymin <- min(es); peak <- if (ymax > -ymin) which.max(es) else which.min(es)
  ggplot(df, aes(rank, es)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    geom_line(colour = col, linewidth = 0.6) +
    geom_segment(aes(x = peak, xend = peak, y = 0, yend = es[peak]), linetype = "dashed",
                 colour = "grey35", linewidth = 0.3) +
    annotate("text", x = peak, y = es[peak], label = sprintf("ES = %.3f", es[peak]),
             hjust = ifelse(es[peak] > 0, -0.05, 1.05), vjust = 0.5, size = 2.8) +
    labs(x = "Ranked gene list", y = "Enrichment score") +
    theme(plot.title = element_blank())
}

## p 值标签：极小值写 < 0.001，避免图上出现 "= 0"
fmt_p_lab <- function(p, prefix = "Log-rank p") {
  if (length(p) == 0 || is.na(p)) return(paste(prefix, "= NA"))
  if (p < 0.001) paste0(prefix, " < 0.001") else paste0(prefix, " = ", formatC(p, format = "f", digits = 3))
}

## KM（删失标记 + log-rank 双侧 p + 每组 n）
km_pub <- function(d, time = "futime", status = "fustat", group = "Risk", xlab = "Time (years)",
                   ylab = "Overall survival probability", pal = c("#0072B2", "#D55E00"),
                   labs = c("Low risk", "High risk")) {
  suppressMessages({library(survival); library(survminer)})
  d[[group]] <- factor(d[[group]], levels = c("low", "high"), labels = labs)
  f <- as.formula(sprintf("Surv(%s, %s) ~ %s", time, status, group))
  s <- do.call(survfit, list(formula = f, data = d))
  sd <- survdiff(f, data = d); pv <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
  g <- ggsurvplot(s, data = d, pval = fmt_p_lab(pv), pval.size = 3.6, pval.coord = c(0.02, 0.12),
                  risk.table = TRUE, risk.table.height = 0.26, censor.size = 2, size = 0.75,
                  palette = pal, legend.title = "", legend.labs = labs,
                  xlab = xlab, ylab = ylab, ggtheme = theme_pub())
  g
}

## 时间依赖 ROC（timeROC）：阶梯曲线 + AUC 95%CI 写进图例
roc_step <- function(d, times = c(1, 3, 5), marker_col = "riskScore",
                     time_col = "futime", status_col = "fustat", n.sim = 500) {
  suppressMessages(library(timeROC))
  r <- timeROC(T = d[[time_col]], delta = d[[status_col]], marker = d[[marker_col]],
               cause = 1, times = times, iid = TRUE)
  ci <- confint(r, parm = NULL, n.sim = n.sim); ci$CI <- ci$CI / 100
  labs <- sprintf("%s-year AUC = %.3f (%.3f-%.3f)", times, r$AUC, ci$CI[, 1], ci$CI[, 2])
  names(labs) <- paste0(times, "-year")
  df <- do.call(rbind, lapply(seq_along(times), function(i) {
    data.frame(fpr = r$FP[, i], tpr = r$TP[, i], Time = names(labs)[i])
  }))
  df$Time <- factor(df$Time, levels = names(labs), labels = labs)
  ggplot(df, aes(fpr, tpr, colour = Time)) +
    geom_step(linewidth = 0.75, direction = "hv") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    scale_colour_manual(values = OI[seq_along(times)], name = "") +
    scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
    labs(x = "1 - Specificity", y = "Sensitivity") + coord_equal() +
    theme(legend.position = c(0.60, 0.20))
}

## 时间依赖 AUC 随时间变化（含 95%CI 带状）
auc_time <- function(d, marker_col = "riskScore", time_col = "futime", status_col = "fustat",
                     tmax = 5, by = 0.25, n.sim = 300) {
  suppressMessages(library(timeROC))
  ts <- seq(by, tmax, by = by)
  r <- timeROC(T = d[[time_col]], delta = d[[status_col]], marker = d[[marker_col]],
               cause = 1, times = ts, iid = TRUE)
  ci <- confint(r, parm = NULL, n.sim = n.sim); ci$CI <- ci$CI / 100
  df <- data.frame(time = ts, auc = as.numeric(r$AUC), lo = ci$CI[, 1], hi = ci$CI[, 2])
  ggplot(df, aes(time, auc)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#0072B2", alpha = 0.16, colour = NA) +
    geom_line(colour = "#0072B2", linewidth = 0.8) +
    geom_point(colour = "#0072B2", size = 1.3) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.5, 1, 0.1)) +
    labs(x = "Time (years)", y = "Time-dependent AUC") + theme(legend.position = "none")
}

## 阶梯状 ROC（含 AUC 95%CI），图例按性能降序
step_roc <- function(roc_list, xlab = "1 - Specificity", ylab = "Sensitivity", digits = 3) {
  if (is.null(names(roc_list)) || any(names(roc_list) == "")) names(roc_list) <- paste0("Model ", seq_along(roc_list))
  o <- order(-vapply(roc_list, function(x) x$auc, 1))
  roc_list <- roc_list[o]
  d <- do.call(rbind, lapply(names(roc_list), function(nm) {
    r <- roc_list[[nm]]
    data.frame(fpr = 1 - r$specificities, tpr = r$sensitivities, Model = nm)
  }))
  d$Model <- factor(d$Model, levels = names(roc_list))
  ggplot(d, aes(fpr, tpr, colour = Model)) +
    geom_step(linewidth = 0.7, direction = "hv") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    scale_colour_manual(values = OI[seq_along(roc_list)], name = "") +
    labs(x = xlab, y = ylab) + coord_equal()
}
