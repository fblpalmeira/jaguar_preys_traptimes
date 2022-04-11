##=========================================
## This program produces Fig 2 in the paper
##=========================================

    ##--------------
    ## Read data etc
    ##--------------
source("jaguar_preys_traptimes_setup.r")

    ##----------------------------
    ## Calculate overlap estimates
    ##----------------------------

cat("Tayassu\n")
ovl2 = ovlest.robust2(times[[1]], times[[5]], kmax=3, c=1)
ovl3 = ovlest.robust2(times[[1]], times[[5]], kmax=3, c=1)
ovl4 = ovlest.robust2(times[[1]], times[[5]], kmax=3, c=1)
print(ovl2)
print(ovl3)
print(ovl4)

# visualizar o tamanho das amostras:
length(times[[1]])
length(times[[5]])
# Calcular as estimativas de sobreposicao:
(Dhats <- overlapEst(times[[1]], times[[5]]))
# Do smoothed bootstrap:
# Generate 999 data sets:
ovl.j <- resample(times[[1]], 999)
ovl.t <- resample(times[[5]], 999)
# Analyse with bootEst, estimating only Dhat1:
bsOut2 <- bootEst(ovl.j, ovl.t) # Takes a while
colMeans(bsOut2)
# Convert column 1 to a vector:
bs2 <- as.vector(bsOut2)
# Get confidence intervals:
bootCI(Dhats[1], bs2)['norm0', ]
bootCI(Dhats[1], bs2)['basic0', ] # escolhido!

cat("Pecari\n")
ovl2 = ovlest.robust2(times[[1]], times[[6]], kmax=3, c=1)
ovl3 = ovlest.robust2(times[[1]], times[[6]], kmax=3, c=1)
ovl4 = ovlest.robust2(times[[1]], times[[6]], kmax=3, c=1)
print(ovl2)
print(ovl3)
print(ovl4)

# visualizar o tamanho das amostras:
length(times[[1]])
length(times[[6]])
# Calcular as estimativas de sobreposicao:
(Dhats <- overlapEst(times[[1]], times[[6]]))
# Do smoothed bootstrap:
# Generate 999 data sets:
ovl.j <- resample(times[[1]], 999)
ovl.pe <- resample(times[[6]], 999)
# Analyse with bootEst, estimating only Dhat1:
bsOut2 <- bootEst(ovl.j, ovl.pe) # Takes a while
colMeans(bsOut2)
# Convert column 1 to a vector:
bs2 <- as.vector(bsOut2)
# Get confidence intervals:
bootCI(Dhats[1], bs2)['norm0', ]
bootCI(Dhats[1], bs2)['basic0', ] # escolhido!

cat("Tapir\n")
ovl2 = ovlest.robust2(times[[1]], times[[7]], kmax=3, c=1)
ovl3 = ovlest.robust2(times[[1]], times[[7]], kmax=3, c=1)
ovl4 = ovlest.robust2(times[[1]], times[[7]], kmax=3, c=1)
print(ovl2)
print(ovl3)
print(ovl4)

# visualizar o tamanho das amostras:
length(times[[1]])
length(times[[7]])
# Calcular as estimativas de sobreposicao:
(Dhats <- overlapEst(times[[1]], times[[7]]))
# Do smoothed bootstrap:
# Generate 999 data sets:
ovl.j <- resample(times[[1]], 999)
ovl.a <- resample(times[[7]], 999)
# Analyse with bootEst, estimating only Dhat1:
bsOut2 <- bootEst(ovl.j, ovl.a) # Takes a while
colMeans(bsOut2)
# Convert column 1 to a vector:
bs2 <- as.vector(bsOut2)
# Get confidence intervals:
bootCI(Dhats[1], bs2)['norm0', ]
bootCI(Dhats[1], bs2)['basic0', ] # escolhido!

cat("Mazama\n")
ovl2 = ovlest.robust2(times[[1]], times[[8]], kmax=3, c=1)
ovl3 = ovlest.robust2(times[[1]], times[[8]], kmax=3, c=1)
ovl4 = ovlest.robust2(times[[1]], times[[8]], kmax=3, c=1)
print(ovl2)
print(ovl3)
print(ovl4)

# visualizar o tamanho das amostras:
length(times[[1]])
length(times[[8]])
# Calcular as estimativas de sobreposicao:
(Dhats <- overlapEst(times[[1]], times[[8]]))
# Do smoothed bootstrap:
# Generate 999 data sets:
ovl.j <- resample(times[[1]], 999)
ovl.m <- resample(times[[8]], 999)
# Analyse with bootEst, estimating only Dhat1:
bsOut2 <- bootEst(ovl.j, ovl.m) # Takes a while
colMeans(bsOut2)
# Convert column 1 to a vector:
bs2 <- as.vector(bsOut2)
# Get confidence intervals:
bootCI(Dhats[1], bs2)['norm0', ]
bootCI(Dhats[1], bs2)['basic0', ] # escolhido!

    ##-----------
    ## Draw Fig 
    ##-----------

png(file = "jaguar_preys_traptimes.png", width = 800, height = 700)

kmax = 3
par(mfrow=c(2,2), mar=c(7,5,5,5), oma=c(1.5,0.5,0.1,0.1), cex.axis=1.5, cex.lab=2)

    cval = 1

    ztig = kplot(times[[1]], kmax, cval)

    z = kplot(times[[5]], kmax, cval)
    plot(z$x, z$y/24, type="l", ylim=c(0,0.16), 
         xlab= "Time of day", ylab="Density of activity", xaxt="n")
    rug(as.vector(times[[5]])/(2*pi))
    zmin = pmin(z$y, ztig$y)/24
    polygon(c(0,ztig$x,1,0), c(0,zmin,0,0), col="lightblue", density=-1, border=NA)
    lines(z$x, z$y/24)
    lines(ztig$x, ztig$y/24, lty="dashed", col="blue")
    axis(1, at=c(0,0.25,0.5,0.75,1), labels=c("0:00", "6:00", "12:00", "18:00", "23:59"))
    text(0.02, 0.15, expression(hat(Delta)[1]), adj=c(0,0))
    text(0.084, 0.152, "= 0.62 (0.50 - 0.76)", adj=c(0,0))
    legend("topright", legend=c("White-lipped peccary", "Jaguar"),
           col=c("black", "blue"), lty=1:2, cex=1.1)
    
    z = kplot(times[[6]], kmax, cval)
    plot(z$x, z$y/24, type="l", ylim=c(0,.16), 
         xlab= "Time of day", ylab="Density of activity", xaxt="n")
    rug(as.vector(times[[6]])/(2*pi))
    zmin = pmin(z$y, ztig$y)/24
    polygon(c(0,ztig$x,1,0), c(0,zmin,0,0), col="lightblue", density=-1, border=NA)
    lines(z$x, z$y/24)
    lines(ztig$x, ztig$y/24, lty="dashed", col="blue")
    axis(1, at=c(0,0.25,0.5,0.75,1), labels=c("0:00", "6:00", "12:00", "18:00", "23:59"))
    text(0.02, 0.15, expression(hat(Delta)[1]), adj=c(0,0))
    text(0.084, 0.152, "= 0.55 (0.40 - 0.65)", adj=c(0,0))
    legend("topright", legend=c("Collared peccary", "Jaguar"),
           col=c("black", "blue"), lty=1:2, cex=1.1)
    
    z = kplot(times[[7]], kmax, cval)
    plot(z$x, z$y/24, type="l", ylim=c(0,0.16), 
         xlab= "Time of day", ylab="Density of activity", xaxt="n")
    rug(as.vector(times[[7]])/(2*pi))
    zmin = pmin(z$y, ztig$y)/24
    polygon(c(0,ztig$x,1,0), c(0,zmin,0,0), col="lightblue", density=-1, border=NA)
    lines(z$x, z$y/24)
    lines(ztig$x, ztig$y/24, lty="dashed", col="blue")
    axis(1, at=c(0,0.25,0.5,0.75,1), labels=c("0:00", "6:00", "12:00", "18:00", "23:59"))
    text(0.02, 0.15, expression(hat(Delta)[1]), adj=c(0,0))
    text(0.084, 0.152, "= 0.54 (0.37 - 0.65)", adj=c(0,0))
    legend("topright", legend=c("Lowland tapir", "Jaguar"),
           col=c("black", "blue"), lty=1:2, cex=1.1)
    
    z = kplot(times[[8]], kmax, cval)
    plot(z$x, z$y/24, type="l", ylim=c(0,0.16), 
         xlab= "Time of day", ylab="Density of activity", xaxt="n")
    rug(as.vector(times[[8]])/(2*pi))
    zmin = pmin(z$y, ztig$y)/24
    polygon(c(0,ztig$x,1,0), c(0,zmin,0,0), col="lightblue", density=-1, border=NA)
    lines(z$x, z$y/24)
    lines(ztig$x, ztig$y/24, lty="dashed", col="blue")
    axis(1, at=c(0,0.25,0.5,0.75,1), labels=c("0:00", "6:00", "12:00", "18:00", "23:59"))
    text(0.02, 0.15, expression(hat(Delta)[1]), adj=c(0,0))
    text(0.084, 0.152, "= 0.68 (0.50 - 0.80)", adj=c(0,0))
    legend("topright", legend=c("Red brocket", "Jaguar"),
           col=c("black", "blue"), lty=1:2, cex=1.1)
 
    
dev.off()

    library(magick)
    library(magrittr) 
    
    # Call back the plot
    plot <- image_read("jaguar_preys_traptimes.png")
    plot2<-image_annotate(plot, "Temporal overlap between jaguar and large preys", 
                          color = "blue", size = 25,
                          location = "10+50", gravity = "north")
    plot3<-image_annotate(plot2, "Data: Palmeira and Trinca (unpublished data) | Visualization by @fblpalmeira 
                          Image credit: T. pecari (Palmeira, FBL), M. americana (Palomo-Munoz, G), Other species (Public Domain) @PhyloPic", 
                          color = "gray", size = 15, 
                          location = "10+50", gravity = "southeast")
    # And bring in a logo
    peccary <- image_read("http://www.phylopic.org/assets/images/submissions/44fb7d4f-6d59-432b-9583-a87490259789.512.png") 
    out1<-image_composite(plot3,image_scale(peccary,"x40"), gravity="north", offset = "+280+110")
    
    tayassu <- image_read("https://raw.githubusercontent.com/fblpalmeira/peccary_abundance/main/img/peccary_avatar.png") 
    tayassu2<-image_flop(tayassu)
    out2<-image_composite(out1,image_scale(tayassu2,"x50"), gravity="north", offset = "-120+110")
 
    mazama <- image_read("http://www.phylopic.org/assets/images/submissions/b5f40112-0cb8-4994-aa70-28ac97ccb83f.512.png") 
    out3<-image_composite(out2,image_scale(mazama,"x60"), gravity="south", offset = "+290+195")
    
    tapir<- image_read("http://www.phylopic.org/assets/images/submissions/8f6b8802-52f9-4f16-8429-0b86ea4a4aa8.512.png") 
    out4<-image_composite(out3,image_scale(tapir,"x50"), gravity="south", offset = "-155+205")
    
    image_browse(out4)
    
    # And overwrite the plot without a logo
    image_write(out2, "jaguar_prey_traptimes2.png")
    