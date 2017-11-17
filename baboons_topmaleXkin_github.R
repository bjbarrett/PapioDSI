#Joan's Baboon Data for Stats

##BRENDAN STARTS

#check code quick for data construction
R
library(rethinking)
library(Cairo)
library(dplyr)

#data is loaded as d
#make charchter or factors integers for stan computation
d$dyad_index <- as.integer(as.factor(d$dyad))
d$focal_index <- as.integer(as.factor(d$focal))
d$partner_index <- as.integer(as.factor(d$partner))
d$group_index <- as.integer(d$focgrp)

###MODEL WITH MALE ONLY (model mMIG in paper)
mmale <- map2stan(
alist(

newdsi ~ dzagamma2( p, mu , scale ),
logit(p) ~ ap + Zp + ap_dyad + bp_male*sametopmale + bpkid*numinfL6mos + bp_phg*PHG ,
log(mu) ~   am + Zm + am_dyad + bm_male*sametopmale + bmkid*numinfL6mos + bm_phg*PHG ,
Zp ~ ap_actor[focal_index] + ap_actor[partner_index],
Zm ~ am_actor[focal_index] + am_actor[partner_index],

	c(ap_actor,am_actor)[focal_index] ~ dmvnorm2( 0 , sigma_actor , Rho_actor ),
	c(ap_dyad,am_dyad)[dyad_index] ~ dmvnorm2( 0 , sigma_dyad , Rho_dyad ),
	c(Rho_actor,Rho_dyad) ~ dlkjcorr(3),
	c(ap,am,bp_male,bm_male,bpkid,bmkid,bp_phg,bm_phg) ~ dnorm(0,2),
	c(sigma_dyad,sigma_actor) ~ dcauchy(0,2),
	scale ~ dcauchy(0,2)
),

data=d, cores=2 , chains=2 , warmup=1500, iter=3500, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE
)
 pars
write.csv(precis(mmale , depth=2 , digits=2)@output , "precismmale.csv" )

par(mfrow=c(2,1) )
a_foc_z <- matrix(0,1000,length(unique(d$focal_index)))
a_prt_z <- matrix(0,1000,length(unique(d$partner_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))


d.pred <- list(
	focal_index=c(1,1),
	partner_index=c(1,1),
	dyad_index=c(1,1),
	sametopmale=c(0,1),
	#rval=rep(mean(d$rval),2),
	numinfL6mos=rep(mean(d$numinfL6mos),2),
	#agediffstd=mean(d$numinfL6mos),
	momkid=c(0,0),
	sibs=c(0,0),
	aunt=c(0,0),
	distkin=c(0,0),
	gmagch=c(0,0),
	PHG=rep(mean(d$PHG),2)
)


link2 <- link(mmale , n=1000 , data=d.pred, replace=
	list(ap_actor=a_foc_z , ap_dyad=a_dyad_z,am_actor=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred <- (1-link2$p)*link2$mu
median(pred[,1])
HPDI(pred[,1])
median(pred[,2])
HPDI(pred[,2])

cairo_pdf("SameTopMaleEffectsDyad.pdf", height=7,width=7)
par(mar=c(5,.5,.5,.5) )
dens(pred[,1], xlim=c(0,15) , xlab="Dyadic Sociality Index" , ylab="" , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE , cex.lab=1.5 , cex.axis=0.5)
ll <- d$newdsi[d$sametopmale==0]
points(ll, rep(-.01,length(ll)), pch=15 , col=col.alpha("orange", alpha=0.1) , cex=0.75 )

shade( density(pred[,1]) , lim= as.vector(HPDI(pred[,1], prob=0.9999)) , col = col.alpha("orange", 0.5))
shade( density(pred[,2]) , lim= as.vector(HPDI(pred[,2], prob=0.9999)) , col = col.alpha("slateblue", 0.5))
ll <- d$newdsi[d$sametopmale==1]


points(ll, rep(-.03,length(ll)), pch=17 , col=col.alpha("slateblue", alpha=0.1) , cex=0.75 )
	axis(1, at = seq(from=0 , to=15, by = 5) ,tck=-0.02 , labels=T )
	axis(1, at = seq(from=0 , to=15, by = 1) ,tck=-0.01 , labels=F)
	#legend("topright", inset=.001, c("different top male" , "same top male"), col=c(col.alpha("orange", 0.5) , col.alpha("slateblue", 0.5) ) , lty=c(1,2) , pch=c(15,17) , bty="n" , cex=1.5)
abline(v=median(pred[,1]) , lty=1)
abline(v=median(pred[,2]) , lty=2)

legend(6.9,1.23, legend = c("", ""),
       col=c(1,1)  , lty=c(1,2),
       lw=1 , cex=1.5, bty="n")

legend(8.8,1.23,, legend = c("Different Top Male", "Same Top Male"),
       col=c(col.alpha("orange", 0.5) , col.alpha("slateblue", 0.5) ) , pch=c(15,17),
       pt.cex=2 , cex=1.5, bty="n")
       #mtext("DSI (dyadic sociality index)" , side=1 , line=2, outer=FALSE , cex=1.5)

dev.off()





d.pred <- list(
 focal_index=1,
 partner_index=1,
 dyad_index=1,
 sametopmale=1,
 #rval=rep(mean(d$rval),2),
 numinfL6mos=mean(d$numinfL6mos),
 #agediffstd=mean(d$numinfL6mos),
 momkid=0,
 sibs=0,
 aunt=0,
 distkin=0,
 gmagch=0,
 PHG=mean(d$PHG)
 )


link2 <- link(m3 , n=1000 , data=d.pred, replace=
	list(ap_actor=a_foc_z , ap_dyad=a_dyad_z,am_actor=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred <- (1-link2$p)*link2$mu

dens(pred, xlim=c(0,14) , show.HPDI=0.89 , main="same topmale" , xlab="dyadic sociality index")
abline(v=median(pred))
ll <- d$newdsi[d$sametopmale==1]
points(ll, rep(0,length(ll)), pch=19 , col=col.alpha("black", alpha=0.15) )


###MODEL WITH SAMETOPMALE and KIN INTERACTION (mKMIG in paper)

mkinXmale <- map2stan(
alist(

newdsi ~ dzagamma2( p, mu , scale ),

logit(p) ~ ap + Zp + ap_dyad + bp_male*sametopmale + (bp_momkid + bp_momkidXmale*sametopmale)*momkid + (bp_sibs + bp_sibsXmale*sametopmale)*sibs 
+ (bp_aunt + bp_auntXmale*sametopmale)*aunt + (bp_gma + bp_gmaXmale*sametopmale)*gmagch + (bp_distkin+ bp_distkinXmale*sametopmale)*distkin + bpkid*numinfL6mos + bp_phg*PHG ,
log(mu) ~   am + Zm + am_dyad + bm_male*sametopmale + (bm_momkid + bm_momkidXmale*sametopmale)*momkid + (bm_sibs + bm_sibsXmale*sametopmale)*sibs 
+ (bm_aunt + bm_auntXmale*sametopmale)*aunt + (bm_gma + bm_gmaXmale*sametopmale)*gmagch + (bm_distkin+ bm_distkinXmale*sametopmale)*distkin + bmkid*numinfL6mos + bm_phg*PHG  ,
Zp ~ ap_actor[focal_index] + ap_actor[partner_index],
Zm ~ am_actor[focal_index] + am_actor[partner_index],

c(ap_actor,am_actor)[focal_index] ~ dmvnorm2( 0 , sigma_actor , Rho_actor ),
c(ap_dyad,am_dyad)[dyad_index] ~ dmvnorm2( 0 , sigma_dyad , Rho_dyad ),
c(Rho_actor,Rho_dyad) ~ dlkjcorr(3),
c(ap,am,bpkid,bmkid,bm_momkid,bm_sibs,bm_aunt,bm_distkin,bp_momkid,bp_sibs,bp_aunt,bp_distkin,bp_gma,bm_gma,bm_phg,bp_phg) ~ dnorm(0,2),
c(bp_male,bm_male,bp_momkidXmale,bm_momkidXmale,bp_sibsXmale,bm_sibsXmale,bp_auntXmale,bm_auntXmale,bp_gmaXmale,bm_gmaXmale,bp_distkinXmale,bm_distkinXmale) ~ dnorm(0,2),
c(sigma_dyad,sigma_actor) ~ dcauchy(0,2),
scale ~ dcauchy(0,2)
),

data=d, cores=2 , chains=2 , warmup=1500, iter=3500, constraints=list(scale="lower=0"), control=list(adapt_delta=0.99) , WAIC=TRUE
)

write.csv(precis(mkinXmale , depth=2 , digits=2)@output , "precismkinXmale.csv" )


d.pred.nk <- list(
	focal_index=c(1,1),
	partner_index=c(1,1),
	dyad_index=c(1,1),
	sametopmale=c(0,1),
	#rval=rep(mean(d$rval),2),
	numinfL6mos=rep(mean(d$numinfL6mos),2),
	#agediffstd=mean(d$numinfL6mos),
	momkid=c(0,0),
	sibs=c(0,0),
	aunt=c(0,0),
	distkin=c(0,0),
	gmagch=c(0,0),
	PHG=rep(mean(d$PHG),2)
)

d.pred.mom <- list(
	focal_index=c(1,1),
	partner_index=c(1,1),
	dyad_index=c(1,1),
	sametopmale=c(0,1),
	#rval=rep(mean(d$rval),2),
	numinfL6mos=rep(mean(d$numinfL6mos),2),
	#agediffstd=mean(d$numinfL6mos),
	momkid=c(1,1),
	sibs=c(0,0),
	aunt=c(0,0),
	distkin=c(0,0),
	gmagch=c(0,0),
	PHG=rep(mean(d$PHG),2)
)

d.pred.sibs <- list(
	focal_index=c(1,1),
	partner_index=c(1,1),
	dyad_index=c(1,1),
	sametopmale=c(0,1),
	#rval=rep(mean(d$rval),2),
	numinfL6mos=rep(mean(d$numinfL6mos),2),
	#agediffstd=mean(d$numinfL6mos),
	momkid=c(0,0),
	sibs=c(1,1),
	aunt=c(0,0),
	distkin=c(0,0),
	gmagch=c(0,0),
	PHG=rep(mean(d$PHG),2)
)

d.pred.gma <- list(
	focal_index=c(1,1),
	partner_index=c(1,1),
	dyad_index=c(1,1),
	sametopmale=c(0,1),
	#rval=rep(mean(d$rval),2),
	numinfL6mos=rep(mean(d$numinfL6mos),2),
	#agediffstd=mean(d$numinfL6mos),
	momkid=c(0,0),
	sibs=c(0,0),
	aunt=c(0,0),
	distkin=c(0,0),
	gmagch=c(1,1),
	PHG=rep(mean(d$PHG),2)
)

d.pred.aunt <- list(
	focal_index=c(1,1),
	partner_index=c(1,1),
	dyad_index=c(1,1),
	sametopmale=c(0,1),
	#rval=rep(mean(d$rval),2),
	numinfL6mos=rep(mean(d$numinfL6mos),2),
	#agediffstd=mean(d$numinfL6mos),
	momkid=c(0,0),
	sibs=c(0,0),
	aunt=c(1,1),
	distkin=c(0,0),
	gmagch=c(0,0),
	PHG=rep(mean(d$PHG),2)
)

d.pred.dist <- list(
	focal_index=c(1,1),
	partner_index=c(1,1),
	dyad_index=c(1,1),
	sametopmale=c(0,1),
	#rval=rep(mean(d$rval),2),
	numinfL6mos=rep(mean(d$numinfL6mos),2),
	#agediffstd=mean(d$numinfL6mos),
	momkid=c(0,0),
	sibs=c(0,0),
	aunt=c(0,0),
	distkin=c(1,1),
	gmagch=c(0,0),
	PHG=rep(mean(d$PHG),2)
)


link2nk <- link(mkinXmale , n=1000 , data=d.pred.nk, replace=
	list(ap_actor=a_foc_z , ap_dyad=a_dyad_z,am_actor=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred.nk <- (1-link2nk$p)*link2nk$mu
median(pred.nk[,1])
HPDI(pred.nk[,1])
median(pred.nk[,2])
HPDI(pred.nk[,2])
link2mom <- link(mkinXmale , n=1000 , data=d.pred.mom, replace=
	list(ap_actor=a_foc_z , ap_dyad=a_dyad_z,am_actor=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred.mom <- (1-link2mom$p)*link2mom$mu
median(pred.mom[,1])
HPDI(pred.mom[,1])
median(pred.mom[,2])
HPDI(pred.mom[,2])
link2sibs <- link(mkinXmale , n=1000 , data=d.pred.sibs, replace=
	list(ap_actor=a_foc_z , ap_dyad=a_dyad_z,am_actor=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred.sibs <- (1-link2sibs$p)*link2sibs$mu
median(pred.sibs[,1])
HPDI(pred.sibs[,1])
median(pred.sibs[,2])
HPDI(pred.sibs[,2])
link2gma <- link(mkinXmale , n=1000 , data=d.pred.gma, replace=
	list(ap_actor=a_foc_z , ap_dyad=a_dyad_z,am_actor=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred.gma <- (1-link2gma$p)*link2gma$mu
median(pred.gma[,1])
HPDI(pred.gma[,1])
median(pred.gma[,2])
HPDI(pred.gma[,2])
link2aunt <- link(mkinXmale , n=1000 , data=d.pred.aunt, replace=
	list(ap_actor=a_foc_z , ap_dyad=a_dyad_z,am_actor=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred.aunt <- (1-link2aunt$p)*link2aunt$mu
median(pred.aunt[,1])
HPDI(pred.aunt[,1])
median(pred.aunt[,2])
HPDI(pred.aunt[,2])
link2dist <- link(mkinXmale , n=1000 , data=d.pred.dist, replace=
	list(ap_actor=a_foc_z , ap_dyad=a_dyad_z,am_actor=a_foc_z , am_dyad=a_dyad_z), WAIC=TRUE)
pred.dist <- (1-link2dist$p)*link2dist$mu
median(pred.dist[,1])
HPDI(pred.dist[,1])
median(pred.dist[,2])
HPDI(pred.dist[,2])

a_foc_z <- matrix(0,1000,length(unique(d$focal_index)))
a_prt_z <- matrix(0,1000,length(unique(d$partner_index)))
a_dyad_z <- matrix(0,1000,length(unique(d$dyad_index)))

dataplots.diff <- cbind(pred.mom[,1],pred.sibs[,1],pred.gma[,1],pred.aunt[,1],pred.dist[,1],pred.nk[,1])
dataplots.same <- cbind(pred.mom[,2],pred.sibs[,2],pred.gma[,2],pred.aunt[,2],pred.dist[,2],pred.nk[,2])


cairo_pdf("KinEffectsDyadXsametopmale.pdf", height=10,width=7)

par(mfrow=c(6,1), cex=1 , mar=c(0,0,0,0) , oma=c(4,0.5,0.5,0.5))
dataplots.diff <- cbind(pred.mom[,1],pred.sibs[,1],pred.gma[,1],pred.aunt[,1],pred.dist[,1],pred.nk[,1])
dataplots.same <- cbind(pred.mom[,2],pred.sibs[,2],pred.gma[,2],pred.aunt[,2],pred.dist[,2],pred.nk[,2])

sublist <- cbind(d$mom , d$sibs , d$gmagch , d$aunt , d$distkin, d$nonkin)
plot.titles <- c("a) Mother/Daughter","b) Sisters","c) Grandmother/Granddaughter","d) Aunt/Niece","e) Distant Kin","f) Non-kin")

for (i in 1:ncol(dataplots.diff)){
dens(dataplots.diff[,i], xlim=c(0,10) , ylim=c(-.035,1.5) , xlab="" , ylab="" , col="white"  , xaxt='n' ,  yaxt='n', zero.line = FALSE , cex.lab=1.5 , cex.axis=0.5 , adj=0.1)
ll <- d$newdsi[sublist[,i]==1 & d$sametopmale==0]
#points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("orange", alpha=0.5) , cex=0.5 )
points(ll, rep(-.018,length(ll)), pch=15 , col=col.alpha("orange", alpha=0.25) , cex=0.4 )

shade( density(dataplots.diff[,i]) , lim= as.vector(HPDI(dataplots.diff[,i], prob=0.99999)) , col = col.alpha("orange", 0.5))
shade( density(dataplots.same[,i]) , lim= as.vector(HPDI(dataplots.same[,i], prob=0.99999)) , col = col.alpha("slateblue", 0.5))
abline(v=median(dataplots.same[,i]) , lty=2)
abline(v=median(dataplots.diff[,i]) , lty=1)

ll <- d$newdsi[sublist[,i]==1 & d$sametopmale==1]
#points(ll, rep(-.025,length(ll)), pch="|" , col=col.alpha("slateblue", alpha=0.5) , cex=0.5 )
points(ll, rep(-.06,length(ll)), pch=17 , col=col.alpha("slateblue", alpha=0.25) , cex=0.4 )
	
	#axis(1, at = seq(from=0 , to=10, by = 1) ,tck=0.01 , labels=F )
	#axis(1, at = seq(from=0 , to=10, by = 1) ,tck=-0.01 , labels=F )
	mtext( plot.titles[i] , side=3 , line=-1, outer=FALSE , cex=1)

if(i==1){

	legend(6.75,1.6, legend = c("", ""),
       col=c(1,1)  , lty=c(1,2),
       lw=1 , cex=1, bty="n")

legend(7.55,1.6,, legend = c("Different Top Male", "Same Top Male"),
       col=c(col.alpha("orange", 0.5) , col.alpha("slateblue", 0.5) ) , pch=c(15,17),
       pt.cex=2 , cex=1, bty="n")

#legend("topright", inset=.001, c("different top male" , "same top male"), fill=c(col.alpha("orange", 0.5) , col.alpha("slateblue", 0.5) ) , lty=c(1,2)  , bty="n" , cex=1)}
#mtext("DSI (dyadic sociality index)" , side=1 , line=2, outer=FALSE , cex=1.5)
}

}
mtext("Dyadic Sociality Index" , side=1 , line=2.5, outer=TRUE , cex=1.5)
#mtext("posterior density" , side=2 , line=2, outer=TRUE , cex=1.5)
axis(1, at = seq(from=0 , to=10, by = 2) ,tck=-0.1 , labels=T )
axis(1, at = seq(from=1 , to=9, by = 2) ,tck=-0.05 , labels=F )
dev.off()
