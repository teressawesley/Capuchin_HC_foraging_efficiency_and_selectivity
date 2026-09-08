d <- read.csv('generated_data/eff_seq_single_proc_s.csv')
cols <- c('man_hands_duration_s','bite_pull_duration_s','bite_shell_duration_s','hit_surface_duration_s','pound_stone_duration_s','roll_scrub_duration_s')
s <- d[which(d$success == 1), cols]
n <- rowSums(!is.na(s) & s != 0)
counts <- tabulate(n, nbins=6)
hands <- !is.na(s$man_hands_duration_s) & s$man_hands_duration_s != 0
with_hands <- sum(n == 2 & hands)
without_hands <- sum(n == 2 & !hands)
cat('Two techniques: with hands =', with_hands, '; without hands =', without_hands, '\n')
png('outputs/success_technique_counts.png', width=1500, height=950, res=180)
par(mar=c(5.5,4.5,3.5,1), family='sans', las=1)
b <- barplot(counts, names.arg=1:6, col='#377D91', border=NA,
             ylim=c(0,47), xlab='Number of techniques occurring',
             ylab='Number of sequences', main=paste0('Successful sequences (success = 1; N = ', sum(counts), ')'))
rect(b[2]-.5, 0, b[2]+.5, with_hands, col='#D48A32', border=NA)
rect(b[2]-.5, with_hands, b[2]+.5, counts[2], col='#677689', border=NA)
text(b[-2], ifelse(counts[-2]>0,counts[-2]-2,1.5), paste0('n=',counts[-2]),
     col=ifelse(counts[-2]>0,'white','#374151'), font=2)
text(b[2], counts[2]+1.6, paste0('n=',counts[2]), font=2, col='#374151')
text(b[2], with_hands/2, sprintf('n=%d\n(%.1f%%)', with_hands,100*with_hands/counts[2]), col='white',font=2,cex=.9)
text(b[2], with_hands+without_hands/2, sprintf('n=%d (%.1f%%)',without_hands,100*without_hands/counts[2]),col='white',font=2,cex=.65)
legend('topright',legend=c('Hand manipulation + another technique','Two techniques without hand manipulation'),
       fill=c('#D48A32','#677689'),border=NA,bty='n',cex=.75)
mtext('Nonzero, non-NA durations across six processing techniques, including roll/scrub.',
      side=1,line=4,cex=.72,col='#4B5563')
dev.off()
