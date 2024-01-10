g_ic(s1 = 0.0,
     s2 = 0.0,
     s1_t = 0.0, 
     s2_t = 0.0,
     sT = 1)


GxIn <- function(s) 3/4 + 1/8 * s^1 + 1/8 * s^2
GxInt <- function(s) 8/3 * (1/8 + s/4)

GxOut <- function(s) 3/4 + 1/4 * s^1
GxOutt <- function(s) 1



n1 <- function(s1, s2, s1t, s2t) GxIn(s1t) * (s1)^0 * GxOut(s2) * (s2t)^0
s <- 0.0
n1(s,s,s,s)

n2 <- function(s1,s2, s1t, s2t) { 
  
  n1(GxIn(s1t) * GxOutt(s2t), 
     GxIn(s2t) * GxOutt(s1t),
     GxInt(s1t) * GxOut(s2t), 
     GxInt(s1t) * GxOut(s2t)
     )
  
}

n3 <- function(s1,s2, s1t, s2t) { 
  
  n2(GxIn(s1t) * GxOutt(s2t), 
     GxIn(s2t) * GxOutt(s1t),
     GxInt(s1t) * GxOut(s2t), 
     GxInt(s1t) * GxOut(s2t)
  )
  
}

n4 <- function(s1,s2, s1t, s2t) { 
  
  n3(GxIn(s1t) * GxOutt(s2t), 
     GxIn(s2t) * GxOutt(s1t),
     GxInt(s1t) * GxOut(s2t), 
     GxInt(s1t) * GxOut(s2t)
  )
  
}

# testing function iteration ----------------------------------------------


s <- 0
n1(s,s,s,s)
n2(s,s,s,s)
n3(s,s,s,s)
n4(s,s,s,s)

# that all works nicely. 

M <- 30
x <- exp(2 * pi * 1i * (0:(M - 1)) / M)
pdf <- Re(fft(n2(s1 = x, s2 = x, s1t = x, s2t = x) |> Conj(), inverse = TRUE))
pdf[pdf < 0] <- 0
pdf <- pdf / sum(pdf)
cas_dist <- tibble(cascade_size = 0:(length(pdf) - 1), prob = pdf) |> pull(prob)

# now try to write a look that will do that.  -----------------------------

citation()
