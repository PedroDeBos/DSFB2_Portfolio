#This file would contain your password. Using dput and as.raw, you can conceal your password as long as they owners don't have acces to this file

pwd_temp<-dput(charToRaw("YourPasswordHere"))
pwd<-as.raw(pwd_temp)
as.raw(pwd_temp)
