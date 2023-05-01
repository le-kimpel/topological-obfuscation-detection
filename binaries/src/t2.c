#include <stdio.h>

int func(int x){
  while(x < 10){
    x+=1;
  }
  return x;
}

int main(void){
  int x = 3;
  puts("ABC");
  puts("CDF");
  x = func(x);
  printf("%d\n", x);
  
}
