
draw_ga <- function() {
  par(mar=c(0,0,0,0), xaxs='i', yaxs='i', family='Helvetica')
  plot.new()
  plot.window(xlim=c(0, 1), ylim=c(0, 1))
  rect(0, 0, 1, 1, col='#F8FAFC', border=NA)
  title_col <- '#12355B'
  accent <- c('#2563A6', '#2F8F83', '#9A5B9A')
  text(0.5, 0.93, 'CSBJ graphical abstract: TNBC-BIR state decomposition', cex=1.35, font=2, col=title_col)
  text(0.5, 0.885, 'Prespecified locked signatures map context-dependent immune-reactive biology', cex=0.82, col='#374151')
  x0 <- c(0.055, 0.365, 0.675)
  x1 <- c(0.325, 0.635, 0.945)
  headers <- c('Locked TNBC-BIR modules', 'Public cohort context map', 'Interpretation boundary')
  body <- list(
    c('Epithelial / plasticity', 'APC / CXCL9 / HLA-II', 'Suppressive myeloid', 'Proliferation controls'),
    c('I-SPY2 / GSE194040: APC axis pCR+', 'GSE25066_RMA: APC null; prolif. pCR+', 'GSE163882: supportive, prolif.-sensitive', 'METABRIC: recurrence biology'),
    c('State decomposition', 'Context-dependent immune biology', 'Transparent negative validation', 'Not a general pCR marker')
  )
  for (i in 1:3) {
    rect(x0[i], 0.18, x1[i], 0.82, col='white', border='#CBD5E1', lwd=1.4)
    rect(x0[i], 0.76, x1[i], 0.82, col=accent[i], border=NA)
    text((x0[i]+x1[i])/2, 0.79, headers[i], col='white', font=2, cex=0.82)
    ys <- seq(0.66, 0.32, length.out=4)
    for (j in 1:4) {
      rect(x0[i]+0.02, ys[j]-0.035, x1[i]-0.02, ys[j]+0.035, col='#EFF6FF', border='#BFDBFE')
      text((x0[i]+x1[i])/2, ys[j], body[[i]][j], cex=0.67, col='#111827')
    }
  }
  arrows(0.335, 0.5, 0.355, 0.5, length=0.08, lwd=2, col='#64748B')
  arrows(0.645, 0.5, 0.665, 0.5, length=0.08, lwd=2, col='#64748B')
  rect(0.055, 0.055, 0.945, 0.125, col='#E0F2FE', border='#7DD3FC')
  text(0.5, 0.092, 'Take-home: reproducible TNBC state mapping with explicit cohort-context limits', cex=0.82, font=2, col='#075985')
}
png('/Users/baiqiu/Desktop/TNBC-BIR/TNBC_BIR_public_only_submission_branch/submission_v1.8_CSBJ/graphical_abstract/graphical_abstract_CSBJ.png', width=2400, height=1350, res=240, type='quartz')
draw_ga()
dev.off()
pdf('/Users/baiqiu/Desktop/TNBC-BIR/TNBC_BIR_public_only_submission_branch/submission_v1.8_CSBJ/graphical_abstract/graphical_abstract_CSBJ.pdf', width=10, height=5.625, useDingbats=FALSE)
draw_ga()
dev.off()
