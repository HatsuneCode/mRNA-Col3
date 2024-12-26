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
fGSEA = function(gene, sig = NULL, scoreType = 'std', minSize = 2, maxSize = 500, type = NULL, species = 'Mus musculus') {
  # species: Homo sapiens / Mus musculus
  suppressMessages(library(fgsea))
  suppressMessages(library(msigdbr))
  if ('KEGG' %in% type) {
    kegg  = msigdbr(species, 'C2', 'KEGG')
    keggn = setNames(lapply(unique(kegg$gs_name), function(i)
      unique(as.character(kegg$gene_symbol)[kegg$gs_name == i])), unique(kegg$gs_name))
    rm(kegg)
    sig = keggn
  }
  if ('GOBP' %in% type) {
    bp   = msigdbr(species, 'C5', 'BP')
    bpn  = setNames(lapply(unique(bp$gs_name), function(i)
      unique(as.character(bp$gene_symbol)[bp$gs_name == i] )), unique(bp$gs_name))
    rm(bp)
    sig = bpn
  }
  if ('GOCC' %in% type) {
    cc   = msigdbr(species, 'C5', 'CC')
    ccn  = setNames(lapply(unique(cc$gs_name), function(i)
      unique(as.character(cc$gene_symbol)[cc$gs_name == i] )), unique(cc$gs_name))
    rm(cc)
    sig = ccn
  }
  if ('GOMF' %in% type) {
    mf   = msigdbr(species, 'C5', 'MF')
    mfn  = setNames(lapply(unique(mf$gs_name), function(i)
      unique(as.character(mf$gene_symbol)[mf$gs_name == i] )), unique(mf$gs_name))
    rm(mf)
    sig = mfn
  }
  set.seed(1)
  ## run gsea analysis by fgsea
  gsea = fgsea(sig, gene, minSize = minSize, maxSize = maxSize, scoreType = scoreType)
  gsea$gene = unlist(lapply(gsea$leadingEdge, function(i) paste(i, collapse = ', ') ))
  data.frame(gsea[, c('pathway', 'NES', 'ES', 'pval', 'padj', 'gene')])
}
GSEA = do.call(rbind, lapply(c('KEGG', 'GOBP', 'GOMF', 'GOCC'), function(i) {
  message('GSEA: ', i)
  gsea = do.call(rbind, lapply(unique(DEG$type), function(g) {
    message(g)
    tmp = DEG[DEG$type == g,]
    gene = sort(setNames(tmp$avg_log2FC, tmp$gene), T)
    df = fGSEA(gene, type = i)
    df$group = g
    df
  }))
  gsea$type = i
  gsea
}))
write.table(GSEA, '4.GSEA.xls', sep = '\t', row.names = F)

## 5. venn ##
df = DEG[DEG$type == 'UVB.Saline vs NoUVB',]
up1 = df$gene[which(df$avg_log2FC >  log2(1.5) & df$p_val < .05)]
dn1 = df$gene[which(df$avg_log2FC < -log2(1.5) & df$p_val < .05)]
df = DEG[DEG$type == 'UVB+hCOL3A1 vs UVB+Saline',]
up2 = df$gene[which(df$avg_log2FC >  log2(1.5) & df$p_val < .05)]
dn2 = df$gene[which(df$avg_log2FC < -log2(1.5) & df$p_val < .05)]
## venn plot
suppressMessages(library(VennDiagram))
venn.diagram(list('UVB-Saline' = up1, 'UVB-COL3LNP' = dn2), filename = 'venn.1.png', fill = c('red', 'blue'),
             disable.logging = T, height = 1600, width = 1600, margin = .15, category.names = c(' ', ' '))
venn.diagram(list('UVB-Saline' = dn1, 'UVB-COL3LNP' = up2), filename = 'venn.2.png', fill = c('blue', 'red'),
             disable.logging = T, height = 1600, width = 1600, margin = .15, category.names = c(' ', ' '))
