## ==parameter==
step =0.05
freq_t1 <- seq(0, 1,by= step)

freq_p1 = 0.1
freq_p2 = 1- freq_p1

a = 0.9

## ==main==
generation = function(freq_t1,freq_t2,freq_p1,freq_p2){
  preference_matrix = matrix(c(
    1,1,1+a,1,
    1,1,1,1+a,
    1,1,1+a,1,
    1,1,1,1+a
  ),nrow = 4)
  
  frequencies = c(freq_p1*freq_t1,
                  freq_p1*freq_t2,
                  freq_p2*freq_t1,
                  freq_p2*freq_t2
  )
  
  weighted = preference_matrix*frequencies
  standard_weighted = weighted/sum(weighted) #每一代都重新校準
  
  new_freq_t1 = standard_weighted[1]+standard_weighted[3]
  new_freq_t2 = standard_weighted[2]+standard_weighted[4]
  
  return(c(new_freq_t1,new_freq_t2))
}

# blank vector field 
plot(NA, xlim = c(0,1), ylim = c(0,1), main = "Vector Field",xlab = 't1', ylab = 't2')


dt1_list = numeric(length(freq_t1))
dt2_list = numeric(length(freq_t1))
i =1

# 以下使用gpt協助
# 繪製向量場
for (t1 in freq_t1) {
  t2 <- 1 - t1 # t1 + t2 = 1
  
  next_gen <- generation(t1, t2, freq_p1, freq_p2)
  dt1 <- next_gen[1] - t1 
  dt2 <- next_gen[2] - t2
  
  dt1_list[i] = dt1
  dt2_list[i] = dt2
  i = i+1
  
  # 畫出每個點的方向箭頭
  arrows(t1, t2, t1 + dt1 * step, t2 + dt2 * step, length = 0.1)
}


