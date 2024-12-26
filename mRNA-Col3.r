## 1. read counts ##
counts = read.table('raw_counts_All_Samples.txt', sep = '\t', header = T, row.names = 1)
counts = counts[rowSums(counts) > 0,]
RNAseq.checkDupRow = function(expr, method = 'mean') {
  ## process symbols
  gene  = sub('^ENS.*?_', '', rownames(expr))
  dgene = unique(gene[duplicated(gene)])
  ## check duplicated symbols
  if (length(dgene)) {
    expr1 = expr[!gene %in% dgene,]
    expr2 = do.call(rbind, lapply(dgene, function(g) {
      e = expr[gene == g, , drop = F]
      if (method == 'mean')
        t(setNames(data.frame(colMeans(e)), g))
    }))
    expr = rbind(expr1, expr2)
    rm(expr1, expr2)
  }
  rm(gene, dgene)
  ## restore symbol
  rownames(expr) = sub('^ENS.*?_', '', rownames(expr))
  expr
}
counts = RNAseq.checkDupRow(counts)
counts = round(counts)

## 2. DEG ##
group = sub('\\..$', '', colnames(counts))
RNAseq.DESeq2 = function(expr, pos = NULL, neg = NULL, name = NULL, exp_cut = 10) {
  suppressMessages(library(DESeq2))
  ## expr split
  exprP = if (length(pos)) expr[, colnames(expr) %in% pos, drop = F] else 
    expr[, !colnames(expr) %in% neg, drop = F]
  exprN = if (length(neg)) expr[, colnames(expr) %in% neg, drop = F] else 
    expr[, !colnames(expr) %in% pos, drop = F]
  ## expr type
  type = paste(paste(colnames(exprP), collapse = ','), 'vs', paste(colnames(exprN), collapse = ',') )
  if (!length(name)) name = type
  message('DEG: ', type)
  ## condition control ~ treatment
  condition = factor( c(rep('Neg', ncol(exprN)), rep('Pos', ncol(exprP))), c('Neg', 'Pos') )
  ## counts
  expr  = cbind(exprN, exprP)
  expr  = expr[rowSums(expr) > 0, , drop = F]
  if (length(exp_cut))
    expr = expr[apply(expr, 1, function(i) !any(i < exp_cut) ), ]
  ##
  exprP = expr[, condition == 'Pos', drop = F]
  exprN = expr[, condition == 'Neg', drop = F]
  ## meta
  meta = data.frame(row.names = colnames(expr), condition)
  ## DESeq2
  dds = DESeqDataSetFromMatrix(countData = expr, colData = meta, design = ~ condition)
  dds = DESeq(dds)
  dds = data.frame(results(dds), check.names = F)
  ## output
  data.frame(p_val = dds$pvalue, avg_log2FC = dds$log2FoldChange, 
             pct.1 = apply(exprP, 1, function(i) sum(i > 0)/ncol(exprP) ),
             pct.2 = apply(exprN, 1, function(i) sum(i > 0)/ncol(exprN) ),
             p_val_adj = dds$padj, gene = rownames(dds), 
             average = rowMeans(expr), median = apply(expr, 1, median), 
             posAvg = rowMeans(exprP), posMed = apply(exprP, 1, median),
             negAvg = rowMeans(exprN), negMed = apply(exprN, 1, median),
             type = name, 
             upDown = factor(ifelse(dds$log2FoldChange > 0, 'Up', 'Down'), c('Up', 'Down')), 
             row.names = NULL )
}
DEG = do.call(rbind, lapply(c('NoUVB', 'UVB.Saline'), function(n) 
  do.call(rbind, lapply(unique(group), function(p) {
    if (p != n) {
      pos = colnames(counts)[group == p]
      neg = colnames(counts)[group == n]
      RNAseq.DESeq2(counts, pos, neg, paste(p, 'vs', n), exp_cut = 10)
    }
  })) ))
write.table(DEG, '2.DEG.bulk.xls', sep = '\t', row.names = F, quote = F)

## 3. volcano ##
cleanGene = function(x, cleanMT = F, cleanSex = T, value = F) {
  idx = !(grepl('^R[P,p][L,l,S,s][0-9,P,p]', x) | grepl('^ENS[M,G]', x) |
            grepl('^AI[0-9][0-9][0-9]', x) | grepl('^B[B,C][0-9][0-9][0-9]', x) |
            grepl('^Gm[0-9]', x) | grepl('^R[N,n][0-9]', x) | grepl('Rik', x) | grepl('^H[B,b].-', x) |
            grepl('^A[A,C,W][0-9]', x) | grepl('^mt-T|^MT-T', x) | grepl('^LINC[0-9]', x) )
  if (cleanMT)  idx = idx & !grepl('^mt-|^MT-', x)
  if (cleanSex) idx = idx & !x %in% c('Xist', 'XIST', 'Ddx3y', 'DDX3Y', 'Eif2s3y', 'EIF2S3Y', 'Uty', 'UTY', 'Kdm5d', 'KDM5D')
  if (value) x[idx] else idx
}
plot_FP = function(df, logFC = 1, padj = .05, label.logFC = Inf, exprAvg = 0, title = 'DEG FC-Padj', adj = F, pmax = 1e-300) {
  suppressMessages(library(ggplot2))
  suppressMessages(library(ggrepel))
  if (!adj) df$p_val_adj = df$p_val
  ## annot
  df$annot = ifelse(df$avg_log2FC > logFC, 'Up', ifelse(df$avg_log2FC < -logFC, 'Down', 'NS'))
  df$annot[df$p_val_adj > padj | is.na(df$p_val_adj)] = 'NS'
  df$annot[df$average < exprAvg] = 'NS'
  df$annot[df$annot == 'Up']   = paste0('Up (', table(df$annot)['Up'], ')')
  df$annot[df$annot == 'Down'] = paste0('Down (', table(df$annot)['Down'], ')')
  df$annot[df$annot == 'NS']   = paste0('NS (', table(df$annot)['NS'], ')')
  df$annot = factor(df$annot, sort(unique(df$annot), T))
  df$lable = abs(df$avg_log2FC) > label.logFC & !grepl('^NS', df$annot)
  ## color
  annot = levels(df$annot)
  idx   = sapply(c('Up', 'NS', 'Down'), function(i) sum(grepl(i, annot)) )
  color = setNames(c(if (idx[1]) 'red', if (idx[2]) 'black', if (idx[3]) 'blue') , annot)
  ## limit
  lim = max(abs(df$avg_log2FC))
  df$p_val_adj[ -log10(df$p_val_adj) > -log10(pmax) ] = pmax
  ## plot
  p = ggplot(df, aes(avg_log2FC, -log10(p_val_adj))) + 
    geom_point(aes(color = annot)) + xlim(c(-lim, lim)) + theme_bw() +
    geom_hline(yintercept = -log10(padj), linetype = 2) +
    geom_vline(xintercept = c(-logFC, logFC), linetype = 2) +
    geom_text_repel(label = ifelse(df$lable, df$gene, NA), family = 'serif') +
    scale_color_manual(values = color) + 
    labs(x = 'Log2 Fold Change', y = if (adj) expression(-log[10](Padj)) else expression(-log[10](Pvalue)), color = NULL, title = title) +
    theme(text = element_text(size = 18), plot.title = element_text(hjust = .5))
  list(plot = p, df = df)
}
DEG = DEG[cleanGene(DEG$gene),]
cut = log2(1.5)
df = DEG[DEG$type == 'UVB.Saline vs NoUVB',]
p  = plot_FP(df, title = 'UVB+Saline vs No UVB', logFC = cut); p$plot
ggsave('3.UVBvsCtrl.volc.pdf', p$plot, w = 7, h = 5.5)
df = DEG[DEG$type == 'UVB.hCOL3A1 vs UVB.Saline',]
p  = plot_FP(df, title = 'UVB+hCOL3A1 vs UVB+Saline', logFC = cut); p$plot
ggsave('3.UVBCOL3vsUVB.volc.pdf', p$plot, w = 7, h = 5.5)

## 4. GSEA ##
GSEA = do.call(rbind, lapply(c('KEGG', 'GOBP', 'GOMF', 'GOCC'), function(i) {
  message('GSEA: ', i)
  gsea = do.call(rbind, lapply(unique(DEG$type), function(g) {
    message(g)
    tmp = DEG[DEG$type == g,]
    gene = sort(setNames(tmp$avg_log2FC, tmp$gene), T)
    df = fGSEA(gene, type = i, species = 'Mus musculus')
    df$group = g
    df
  }))
  gsea$type = i
  gsea
}))
write.table(GSEA, '4.GSEA.xls', sep = '\t', row.names = F)
