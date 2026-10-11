options(stringsAsFactors = FALSE, warn = -1)
ncase <- 544; nctrl <- 378749
neff <- 4 * ncase * nctrl / (ncase + nctrl)
rows <- list()
for (r2 in c(0.01, 0.03, 0.05, 0.08)) {
  f <- function(b, a) pchisq(qchisq(1 - a, 1), 1, ncp = neff * r2 * b^2, lower.tail = FALSE) - 0.80
  or_bonf <- tryCatch(exp(uniroot(f, c(0.001, 2), a = 0.05 / 49)$root), error = function(e) NA)
  or_005 <- tryCatch(exp(uniroot(f, c(0.001, 2), a = 0.05)$root), error = function(e) NA)
  rows[[length(rows) + 1]] <- data.frame(R2 = r2, min_OR_bonferroni = round(or_bonf, 2), min_OR_p005 = round(or_005, 2))
}
res <- do.call(rbind, rows)
print(res, row.names = FALSE)
write.table(res, "F:/0000_Bioinfo/222_tcgaSHAP_胶质瘤/04_results/audit/mr_power_min_or.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
