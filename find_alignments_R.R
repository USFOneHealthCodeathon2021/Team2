find_alignments_R <- function(seq1, seq2) {
  
  
  # INPUTS:
  # seq1: the original Wuhan sequence as a character string
  # seq2: user input sequence as a character string
  
  # OUTPUTS:
  # a visual graphic showing identical sequence regions and mutations between covid sequences
  
  # PARAMETERS: (hardcoded as these are ideal for Covid sequences):
  min_seq_len = 10;
  match_len = 5;
  indel_len = 0;
  
  
  # convert seq1 and seq2 to uppercase in case they're not already:
  seq1 <- toupper(seq1);
  seq2 <- toupper(seq2);
  
  # get lengths of each sequence:
  len1 <- nchar(seq1);
  len2 <- nchar(seq2);
  
  # find alignment for sequences in both orders:
  reps1 <- make_align(seq1, seq2, min_seq_len, match_len, indel_len);
  reps2 <- make_align(seq2, seq1, min_seq_len, match_len, indel_len);
  
  # combine reps1 and reps2 into one variable, switching the columns in reps2:
  temp <- reps2[,1];
  reps2[,1] <- reps2[,3];
  reps2[,3] <- temp;
  temp <- reps2[,2];
  reps2[,2] <- reps2[,4];
  reps2[,4] <- temp;
  reps <- rbind(reps1,reps2);
  
  # merge the acceptable overlaps in reps:
  num_reps = nrow(reps);
  r1 <- 1;
  while(r1 <= num_reps-1) {
    r2 <- r1 + 1;
    while(r2 <= num_reps) {
      if(reps[r1,1]==reps[r2,1] && reps[r1,2]==reps[r2,2] && reps[r1,3]==reps[r2,3] && reps[r1,4]== reps[r2,4]) {
        reps <- reps[-c(r2:r2),];
        num_reps <- num_reps - 1;
      } else if(reps[r2,1] >= reps[r1,1] && reps[r2,1] <= reps[r1,2] || reps[r2,2] >= reps[r1,1] && reps[r2,2] <= reps[r1,2]) {
        if(reps[r1,1]-reps[r2,1] == reps[r1,3]-reps[r2,3] && reps[r1,2]-reps[r2,2] == reps[r1,4]-reps[r2,4]) {
          reps[r1,1] <- min(reps[r1,1],reps[r2,1]);
          reps[r1,2] <- max(reps[r1,2],reps[r2,2]);
          reps[r1,3] <- min(reps[r1,3],reps[r2,3]);
          reps[r1,4] <- min(reps[r1,4],reps[r2,4]);
          reps <- reps[-c(r2:r2),];
          num_reps <- num_reps - 1;
        }
      } else {
        r2 <- r2 + 1;
      }
    }
    r1 <- r1 + 1;
  }
  
  # recreate reps:
  reps <- matrix(reps,ncol = 4);

  # generate coordinates of identical sequence region boxes to be plotted:
  spacer <- 500;
  spaced_reps <- matrix(reps,ncol = 4);
  num_rows <- nrow(spaced_reps);
  for(r in 1:num_rows) {
    spaced_reps[r,] <- spaced_reps[r,] + spacer*r;
  }
  
  # make a visual display of the sequence:
  plot(0,0,
       main = 'Comparison of input CoVid sequence and original human Wuhan strain',
       xlim = c(0,max(spaced_reps[num_rows,2],spaced_reps[num_rows,4])+spacer),
       ylim = c(0,3),
       xlab = 'base pairs',
       ylab = 'strains',
       xaxt = 'n',
       yaxt = 'n'
       )
  par(new = TRUE);
  
  # add labels for the y-axis:
  axis(side = 2, at = c(0.5,2.5), labels = c('original Wuhan strain','input sequence'));
  
  # plot the other elements:
  for(r in 1:num_rows) {
    # plot rectangles of identical sequence regions:
    rect(spaced_reps[r,1], 0, spaced_reps[r,2], 1, col = rgb(0,0.9,0.3,0.5))
    rect(spaced_reps[r,3], 2, spaced_reps[r,4], 3, col = rgb(0,0.9,0.3,0.5))
    # plot a dotted line between identical regions:
    midpt1 <- spaced_reps[r,1] + (spaced_reps[r,2]-spaced_reps[r,1])/2;
    midpt2 <- spaced_reps[r,3] + (spaced_reps[r,4]-spaced_reps[r,3])/2;
    lines(c(midpt1,midpt2),c(1,2),col = rgb(0,0.4,0),lty = 'dashed')
    # label identical regions with base pair length:
    text(midpt1,0.5,paste(spaced_reps[r,2]-spaced_reps[r,1]+1,'bp'),srt = 90)
    text(midpt1,2.5,paste(spaced_reps[r,4]-spaced_reps[r,3]+1,'bp'),srt = 90)
    # plot mutations:
    # mutations that occur before identical regions:
    if(r==1) {
      if(reps[1,1] > 1) {
        misalign1 <- substr(seq1,1,reps[1,1]-1);
        if(reps[1,1]-1 == 1) {
          bp <- '(1)';
        } else {
          bp <- paste('(1:',reps[1,1]-1,')',sep = '');
        }
        text(spacer/2,0.5,paste(misalign1,bp),srt = 90)
      }
      if(reps[1,3] > 1) {
        misalign2 <- substr(seq2,1,reps[1,3]-1);
        if(reps[1,3]-1 == 1) {
          bp <- '(1)';
        } else {
          bp <- paste('(1:',reps[1,3]-1,')',sep = '');
        }
        text(spacer/2,2.5,paste(misalign2,bp),srt = 90)
      }
    }
    # mutations that occur between identical regions:
    if(r < nrow(reps)) {
      misalign1 <- substr(seq1,reps[r,2]+1,reps[r+1,1]-1)
      if(reps[r+1,1]-reps[r,2] == 2) {
        bp <- paste('(',reps[r,2]+1,')',sep = '');
      } else {
        bp <- paste('(',reps[r,2]+1,':',reps[r+1,1]-1,')',sep = '');
      }
      text(spaced_reps[r,2]+1+(spaced_reps[r+1,1]-1-(spaced_reps[r,2]+1))/2,0.5,paste(misalign1,bp),srt = 90)
      misalign2 <- substr(seq2,reps[r,4]+1,reps[r+1,3]-1)
      if(reps[r+1,3]-reps[r,4] == 2) {
        bp <- paste('(',reps[r,4]+1,')',sep = '');
      } else {
        bp <- paste('(',reps[r,4]+1,':',reps[r+1,3]-1,')',sep = '');
      }
      text(spaced_reps[r,4]+1+(spaced_reps[r+1,3]-1-(spaced_reps[r,4]+1))/2,2.5,paste(misalign2,bp),srt = 90)
    } else if(r == nrow(reps)) {
      if(reps[nrow(reps),2] < len1) {
        misalign1 <- substr(seq1,reps[nrow(reps),2]+1,len1);
        if(reps[r,2] == len1-1) {
          bp <- paste('(',len1,')',sep = '');
        } else {
          bp <- paste('(',reps[r,2]+1,':',len1,')',sep = '');
        }
        text(spaced_reps[nrow(reps),2]+spacer/2,0.5,paste(misalign1,bp),srt = 90)
      }
      if(reps[nrow(reps),4] < len2) {
        misalign2 <- substr(seq2,reps[nrow(reps),4]+1,len2);
        if(reps[r,4] == len2-1) {
          bp <- paste('(',len2,')',sep = '');
        } else {
          bp <- paste('(',reps[r,4]+1,':',len2,')',sep = '');
        }
        text(spaced_reps[nrow(reps),4]+spacer/2,2.5,paste(misalign2,bp),srt = 90)
      }
    }
    
    # keep current figure:
    par(new = TRUE);
  }
}


#######################################################################################################################

make_align <- function(seq1, seq2, min_seq_len, match_len, indel_len) {

  
  # determine the length of each input sequence:
  len1 <- nchar(seq1);
  len2 <- nchar(seq2);
  
  # initialize the variables that will hold alignments for each sequence:
  align1 <- matrix(nrow = 0,ncol = 2);
  align2 <- matrix(nrow = 0,ncol = 2);
  
  # find the ranges where sequences match exactly:
  OS1 <- 0;
  while(OS1 <= len1-match_len) {
    OS2 <- 0;
    while(OS2 <= len2-match_len && OS1 <= len1-match_len) {
      S1 <- substr(seq1,1+OS1,OS1+match_len);
      S2 <- substr(seq2,1+OS2,OS2+match_len);
      if(S1==S2) {
        expanding = TRUE;
        expand = 0;
        while(expanding) {
          if(OS1+match_len+expand+1 <= len1 && OS2+match_len+expand+1 <= len2 && substr(seq1,OS1+match_len+expand+1,OS1+match_len+expand+1)==substr(seq2,OS2+match_len+expand+1,OS2+match_len+expand+1)) {
            expand <- expand+1;
          } else {
            expanding <- FALSE;
          }
        }
        align1 <- rbind(align1,c(1+OS1,OS1+match_len+expand));
        align2 <- rbind(align2,c(1+OS2,OS2+match_len+expand));
        OS1 <- OS1+match_len+expand+1;
        OS2 <- OS2+match_len+expand+1;
      } else {
        OS2 <- OS2 + 1;
      }
    }
    OS1 <- OS1 + 1;
  }
  
  # make a matrix of possible transitions between matching ranges:
  aligns = nrow(align1);
  mat = matrix(nrow = aligns, ncol = aligns);
  for(r in seq(1,aligns,1)) {
    for(c in seq(1,aligns,1)) {
      diff1 <- align1[c,1]-align1[r,2]-1;
      diff2 <- align2[c,1]-align2[r,2]-1;
      if(abs(diff1) <= indel_len && abs(diff2) <= indel_len && abs(diff1-diff2) <= indel_len) {
        mat[r,c] = 1;
      } else {
        mat[r,c] = 0;
      }
    }
  }
  
  # find all possible starting ranges for longer sequences, including lone ranges:
  starts <- matrix(nrow = 1,ncol = 0);
  for(a in 1:aligns) {
    row_sum <- sum(mat[a,]);
    col_sum <- sum(mat[,a]);
    if(row_sum > 0 && col_sum == 0 || row_sum == 0 && col_sum == 0) {
      starts <- c(starts,a);
    }
  }
  
  # make a list of all possible sequences of ranges:
  listo <- list();
  if(length(starts) > 0) {
    for(s in 1:length(starts)) {
      listo[[s]] <- starts[s];
    }
    k <- 1;
    while(k <= length(listo)) {
      len_seq <- length(listo[[k]]);
      last_el = listo[[k]][len_seq];
      tails = c();
      for(c in 1:ncol(mat)) {
        if(mat[last_el,c] != 0) {
          tails <- c(tails,c);
        }
      }
      tail <- 1;
      while(tail <= length(tails)) {
        if(tails[tail] < last_el) {
          tails[tail] <- c();
        } else {
          tail <- tail + 1;
        }
      }
      if(length(tails) > 0) {
        listo[[k]] <- c(listo[[k]],tails[1]);
        if(length(tails) > 1) {
          for(t in 2:length(tails)) {
            #listo(size(listo,2)+1).seq = [listo(k).seq(1:len_seq), tails(t)];
            listo[[length(listo)+1]] <- c(listo[[k]][1:len_seq], tails[t]);
          }
        }
      } else {
        k <- k + 1;
      }
    }
  }
  
  reps <- matrix(nrow = 0, ncol = 4);
  i <- 0;
  if(length(listo) > 0) {
    for(k in 1:length(listo)) {
      first_el <- listo[[k]][1];
      last_el <- listo[[k]][length(listo[[k]])];
      range1 <- c(align1[first_el,1],align1[last_el,2]);
      range2 <- c(align2[first_el,1],align2[last_el,2]);
      if(range1[2]-range1[1]+1 >= min_seq_len && range2[2]-range2[1]+1 >= min_seq_len) {
        i <- i + 1;
        reps <- rbind(reps,c(range1[1],range1[2],range2[1],range2[2]));
      }
    }
  }
  
  return(reps)
  
}





