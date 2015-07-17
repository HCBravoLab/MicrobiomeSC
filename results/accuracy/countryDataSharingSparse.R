alpha = .05
m1 = read.csv("sparse_GenDAKenya")
r1 = p.adjust(m1[,3], method="fdr")
r1Above = r1[which(r1<=alpha)]

m2 = read.csv("sparse_GenDAGambia")
r2 = p.adjust(m2[,3], method="fdr")
r2Above = r2[which(r2<=alpha)]

m3 = read.csv("sparse_GenDAMali")
r3 = p.adjust(m3[,3], method="fdr")
r3Above = r3[which(r3<=alpha)]

m4 = read.csv("sparse_GenDABangladesh")
r4 = p.adjust(m4[,3], method="fdr")
r4Above = r4[which(r4<=alpha)]

m5 = read.csv("sparse_GenDAKenyaGambia")
r5 = p.adjust(m5[,3], method="fdr")
r5Above = r5[which(r5<=alpha)]

m6 = read.csv("sparse_GenDAKenyaMali")
r6 = p.adjust(m6[,3], method="fdr")
r6Above = r6[which(r6<=alpha)]

m7 = read.csv("sparse_GenDAKenyaBangladesh")
r7 = p.adjust(m7[,3], method="fdr")
r7Above = r7[which(r7<=alpha)]

m8 = read.csv("sparse_GenDAMaliGambia")
r8 = p.adjust(m8[,3], method="fdr")
r8Above = r8[which(r8<=alpha)]

m9 = read.csv("sparse_GenDAGambiaBangladesh")
r9 = p.adjust(m9[,3], method="fdr")
r9Above = r9[which(r9<=alpha)]

m10 = read.csv("sparse_GenDAMaliBangladesh")
r10 = p.adjust(m10[,3], method="fdr")
r10Above = r10[which(r10<=alpha)]

length(r1Above)
length(r2Above)
length(r3Above)
length(r4Above)
length(r5Above)
length(r6Above)
length(r7Above)
length(r8Above)
length(r9Above)
length(r10Above)

