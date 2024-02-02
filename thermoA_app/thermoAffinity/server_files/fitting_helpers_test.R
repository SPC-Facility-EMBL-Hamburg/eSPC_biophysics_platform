rm(list=ls())

scripts_path <- "/home/osva/mst_shiny/server_files/fitting_helpers/"
scripts <- list.files(scripts_path)

lapply(paste0(scripts_path,scripts), source)

lig_test <- 3000/(2**(16:6))

frac_list <- lapply(lig_test, function(x)frac_general_two_sites(5,5,x,1,10))

f_list    <- sapply(frac_list, function(x) 203*(x$ab + x$aba + x$ba) + 120*x$b) 
f_list    <- f_list + rnorm(length(f_list),0,30)

df_test <- data.frame("x"=lig_test,y=f_list)

ggplot(df_test,aes(x,y))+
  geom_point()+scale_x_continuous(trans = "log")
  
rm(fit)

fit <- fit_two_sites_one_kd_shared_signal(df_test$y,df_test$x,10)


fit$tidy_fit

y_hat <- predict(fit$fit_obj,df_test$x)

df_test_pred <- data.frame("x"=lig_test,"y"=y_hat)

ggplot(df_test,aes(x,y))+
  geom_point()+
  scale_x_continuous(trans = "log")+
  geom_point(data=df_test_pred,aes(x,y),color="red")
 
#  <td></td>
#  <td>Two Sites, One Kd, Different Signal</td>
#  <td>Two Sites, One Kd, Shared Signal, Cooperativity</td>
#  <td>Two Sites, One Kd, Different Signal, Cooperativity</td>
#  <td>Two Sites, Two Kds, Shared Signal</td>
#  <td>Two Sites, Two Kds, Different Signal</td>
