/*****************************************
[補足]
*****************************************/

#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <sys/stat.h>
#include <algorithm>

using vd = std::vector<double>;
using vvd = std::vector<vd>;
using vi = std::vector<int>;
using vvi = std::vector<vi>;

/**********************config*************************/
// 0:熱浴法 1:メトロポリス法
const int METHOD = 0;
// 0:ランダム 1:格子状 2:網状 3:基底状態
const int INIT = 0;
/*****************************************************/

/******************計算条件********************/
// 要素数
const int n = 10000;

/**********************一様乱数の定義*************************/
//非決定的な乱数生成機、これを使ってシードをランダムに選ぶ
std::random_device rnddev;
//シード値の指定に使う乱数
std::mt19937 mt(rnddev());
// 一様乱数
std::uniform_real_distribution<double> rndd(0.0,1.0);
// 平均 2.5, 標準偏差 1.0 の正規分布
std::normal_distribution<> rndn(2.5, 1.0);
/************************************************************/

void jackknife(vd &A);

int main(){
  vd A_uni(n), A_normal(n);
  for(int i = 0; i < n; i++) {
    A_uni[i] = rndd(mt);
    A_normal[i] = rndn(mt);
  }

}

void jackknife(vd &A){
  int M = (int)A.size();
  int M_m; // bin の数
  double jsample1, jsample2, sum1, sum2;
  vd S(M+1,0), S2(M+1,0), bin_size;
  // 累積和を計算
  for(int i = 0; i < M; i++){
    S[i+1] = S[i] + A[i];
    S2[i+1] = S2[i] + A[i]*A[i];
  }
  // M の約数を列挙して bin サイズを決定
  for(int i = 1; i*i <= n; i++) {
    if(M%i == 0){
      bin_size.push_back(i);
      if(i*i == M) continue;
      bin_size.push_back(M%i);
    }
  }
  // 一応降順にソートしておく
  std::sort(bin_size.begin(), bin_size.end());
  // M の各約数 m について jackknife 方を適用
  for(int m : bin_size) {
    // m 分割された各グループが持つデータの個数が M_m
    M_m = M/m;
    sum1 = 0.0, sum2 = 0.0;
    for(int i = 0; i < M; i++) {
      if(i%M_m == 0){
        jsample1 = 0.0;
        jsample2 = 0.0;
      }
      jsample1 += v[i];
      jsample2 += v[i]*v[i];
    }
  }
}