##mix colours

##takes in as named colour only!
mixColours <- function(col1,col2,col3=NULL,alphaVal=NULL){

  if(is.null(alphaVal)){
    alphaVal <- 255
  }

  cv1 <- t(col2rgb(col1))
  cv2 <- t(col2rgb(col2))

  if(! is.null(col3)){
    cv3 <- t(col2rgb(col3))
    cv1 <- (cv1 + cv3)/2
  }

  return(rgb(alpha=255,
  red=round(sum(cv1[1],cv2[1])/2),
  green=round(sum(cv1[2],cv2[2])/2),
  blue=round(sum(cv1[3],cv2[3])/2),maxColorValue=255))
}
