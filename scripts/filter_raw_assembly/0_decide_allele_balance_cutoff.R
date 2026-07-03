#assuming average depth is 34


bts <- cbind.data.frame(AC=1:9, pval=sapply(1:9,function(x){
    bt<-binom.test(x,10,p=0.5,alternative = "two.sided")
    return(bt$p.value)}))


bts$filter1<- bts$pval < 0.001
bts$filter2<- bts$pval < 0.0001
bts$filter3<- bts$pval < 0.00001
bts$filter4<- bts$pval < 0.000001
colSums(bts)

cbind(bts$AC,bts$filter2)
5/34
