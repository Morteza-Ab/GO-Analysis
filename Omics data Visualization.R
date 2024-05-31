install.packages('GOplot')
library(GOplot)
data(EC)
head(EC$david)

circ <- circle_dat(EC$david, EC$genelist)
GOBar(subset(circ, category == 'BP'))
GOBar(circ, display = 'multiple')
GOBar(circ, display = 'multiple', title = 'Z-score coloured barplot', zsc.col = c('yellow', 'black', 'cyan'))

GOBubble(circ, labels = 3)
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), display = 'multiple', labels = 3)  
GOCircle(circ)

chord <- chord_dat(circ, EC$genes, EC$process)
head(chord)
chord <- chord_dat(data = circ, genes = EC$genes)
chord <- chord_dat(data = circ, process = EC$process)
chord <- chord_dat(data = circ, genes = EC$genes, process = EC$process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)




GOCluster(circ, EC$process, clust.by = 'logFC', term.width = 2)


l1 <- subset(circ, term == 'heart development', c(genes,logFC))
l2 <- subset(circ, term == 'plasma membrane', c(genes,logFC))
l3 <- subset(circ, term == 'tissue morphogenesis', c(genes,logFC))
GOVenn(l1,l2,l3, label = c('heart development', 'plasma membrane', 'tissue morphogenesis'))
