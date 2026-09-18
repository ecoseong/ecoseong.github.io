---
layout: post
title: "[투자론 2장] 확률분포에서 주가모형과 통계적 추정까지"
date: 2026-09-18 09:10:00 +0900
categories: [Finance, Investments]
tags: [probability, brownian-motion, volatility, statistics]
math: true
toc: true
description: "정규분포, 위너 과정, 기하 브라운 운동, 두꺼운 꼬리, 표본통계량과 표준오차를 수식으로 전개한다."
---

이 글은 MIT 15.433 Investments의 **Class 2: Securities, Random Walk on Wall Street**(Reto r. Gallati, 2003년 2월 5일)을 바탕으로 한 학습 노트다. 단순한 정규모형으로 수익률을 설명한 뒤, 그 모형이 놓치는 위험과 추정의 불확실성까지 살펴본다.

원문의 수학적 오기는 수정했으며 글 끝에 정오표를 첨부했다. **이토 보조정리, GBM의 정확한 해, 자기상관을 반영한 표준오차 등은 보충 내용**이다. 긴 유도는 이 글 끝의 부록에 정리했다.

**시리즈:** [1장]({% post_url 2026-09-18-investments-01-financial-system %}) · 2장 · [이 장의 부록](#chapter-appendix)

## 원문 기호와 이 글의 약속

| 원문 기호 | 이 글에서의 의미 |
|---|---|
| $$z,\ dz,\ \Delta z,\ z(T)-z(0)$$ | 위너 과정과 그 증분. 원문 기호를 그대로 사용 |
| $$x,\ dx,\ a,\ b$$ | 일반화 위너 과정과 드리프트·확산 계수 |
| $$S,\ S_0,\ S_T,\ \mu,\ \sigma$$ | 주가와 GBM의 드리프트·변동성. $$S$$는 현재 시점의 가격 |
| $$r,\ r_i$$ | 수익률. 단순수익률인지 로그수익률인지 각 절에서 명시 |
| $$X=(r-\mu)/\sigma$$ | 표준화한 수익률 |
| $$\mu,\ \sigma^2$$, 표본통계 절 | 원문처럼 표본평균·경험분산을 뜻함 |
| $$\mathrm{skew},\ \mathrm{kurt}$$ | 왜도·첨도. 기호는 유지하고 분모 오류는 수정 |
| $$N(\mu,\sigma)$$ | 원문 식 (29)의 평균·표준편차 관례 |

원문은 서로 다른 문맥에서 같은 기호를 재사용한다. 이를 그대로 따르되 각 절에서 정의를 분명히 한다. 증명을 위해 모집단과 표본을 함께 비교할 때의 $$\mu_{\mathrm{pop}},\sigma_{\mathrm{pop}}$$, 불편분산 $$s^2$$, 두 종류 수익률을 함께 적을 때의 $$\mathrm{simple},\mathrm{log}$$ 윗첨자 등은 **보충 설명용 표시**다. 원문에 대응하는 기본 식의 기호를 대체하지 않는다.

## 1. 확률분포는 결과와 가능성을 함께 표현한다

미래 결과가 확정되지 않은 변수를 확률변수라고 한다. 확률분포는 가능한 값과 그 값이 발생할 가능성을 함께 나타낸다.

### 1.1 베르누이분포와 이항분포

원문과 같이 $$X\in\{0,1\}$$이고 $$\Pr(X=0)=p$$, $$\Pr(X=1)=1-p$$로 둔다.

$$
\Pr(X=x)=p^{1-x}(1-p)^x,\qquad x\in\{0,1\}.
$$

$$X^2=X$$이므로

$$
\mathbb E[X]=1-p,\qquad
\operatorname{Var}(X)=\mathbb E[X^2]-\mathbb E[X]^2
=(1-p)-(1-p)^2=p(1-p).
$$

보충적으로 독립 시행 $$n$$번에서 **1이 나온 횟수**를 $$Y$$라고 하면

$$
Y=\sum_{i=1}^{n}X_i\sim\operatorname{Binomial}(n,1-p),
\qquad \Pr(Y=k)=\binom nk(1-p)^kp^{n-k}.
$$

원문은 이 두 값의 예를 이항분포라고 소개한다. 더 구체적으로는 베르누이분포이며, 시행 횟수가 1인 이항분포다. 여기서 $$p$$는 **0의 확률**이라는 원문 정의를 유지했다.

### 1.2 연속분포에서 밀도는 확률 자체가 아니다

연속확률변수의 확률밀도함수가 $$f(x)$$라면

$$
\Pr(a\le X\le b)=\int_a^bf(x)\,dx,
\qquad \int_{-\infty}^{\infty}f(x)\,dx=1.
$$

밀도의 높이 $$f(x)$$가 아니라 **구간 아래의 면적**이 확률이다. 한 점의 확률은 0이며, 밀도는 1보다 클 수도 있다.

## 2. 정규분포와 표준화

원문 식 (29)에 맞춰 이 글의 정규분포 표기는 **평균과 표준편차** 순서로 쓴다.

$$
\boxed{r\sim N(\mu,\sigma)\quad\text{에서 두 번째 인수는 표준편차이다.}}
$$

일반 통계 교재의 $$N(\mu,\sigma^2)$$ 표기는 두 번째 인수가 분산인 다른 관례다. 여기서는 원문과 대조하기 쉽게 표준편차 관례를 일관되게 사용한다. 따라서 $$\sigma>0$$일 때

$$
f_r(r)=\frac{1}{\sigma\sqrt{2\pi}}
\exp\left[-\frac{(r-\mu)^2}{2\sigma^2}\right].
$$

원문 식 (30)처럼 표준화한 변수를 $$X$$로 쓰면

$$
X=\frac{r-\mu}{\sigma}\sim N(0,1).
$$

원문의 $$N(x)$$는 표준정규 누적분포함수도 뜻한다. **인수 하나인 $$N(x)$$는 누적확률, 인수 둘인 $$N(\mu,\sigma)$$는 분포 표기**라는 점을 구분한다.

$$
\Pr(r\le c)=N\left(\frac{c-\mu}{\sigma}\right),
$$

$$
\Pr(a\le r\le b)
=N\left(\frac{b-\mu}{\sigma}\right)
-N\left(\frac{a-\mu}{\sigma}\right).
$$

| 평균 주변 구간 | 정규분포의 포함 확률 |
|---|---:|
| $$\mu\pm\sigma$$ | 약 68.27% |
| $$\mu\pm2\sigma$$ | 약 95.45% |
| $$\mu\pm3\sigma$$ | 약 99.73% |

이 숫자는 모든 분포에 공통으로 적용되는 법칙이 아니다.

## 3. 정규분포에서 확률과정으로

정규분포 하나는 특정 시점의 결과를 묘사한다. 주가처럼 시간에 따라 움직이는 대상을 표현하려면 확률변수의 집합인 **확률과정**이 필요하다.

이산시간 랜덤워크의 예는

$$
X_{k+1}=X_k+\varepsilon_{k+1}
$$

이다. 브라운 운동은 연속시간 확률과정으로, 적절히 시간·공간을 조정한 랜덤워크의 극한으로 이해할 수 있다. 이산시간 랜덤워크와 연속시간 브라운 운동은 관련되어 있지만 동일한 정의는 아니다.

원문은 Bachelier의 연구에서 시작하여 브라운 운동이 금융의 가격모형과 옵션이론에 연결된 역사를 소개한다. 여기서는 모형의 수학적 구조에 집중한다.

## 4. 위너 과정: 변동의 크기가 왜 시간의 제곱근인가?

### 4.1 정의

표준 위너 과정 $$z(t)$$는 다음 성질을 갖는다.

1. \$$z(0)=0$$.
2. 겹치지 않는 시간구간의 증분은 독립이다.
3. \$$z(t+h)-z(t)\sim N(0,\sqrt h)$$ for \$$h>0$$.
4. 경로가 거의 확실하게 연속이다.

작은 구간 $$\Delta t$$에서는

$$
\Delta z=\varepsilon\sqrt{\Delta t},\qquad \varepsilon\sim N(0,1)
$$

로 표현할 수 있다. 따라서

$$
\mathbb E[\Delta z]=0,
\quad \operatorname{Var}(\Delta z)=\Delta t,
\quad \operatorname{SD}(\Delta z)=\sqrt{\Delta t}.
$$

### 4.2 긴 기간의 분산 — 짧은 증명

$$T=N\Delta t$$인 구간을 나누면

$$
z(T)-z(0)=\sum_{i=1}^{N}\varepsilon_i\sqrt{\Delta t}.
$$

기대값의 선형성으로

$$
\mathbb E[z(T)-z(0)]=\sum_{i=1}^{N}0=0.
$$

독립성으로 공분산이 0이므로

$$
\begin{aligned}
\operatorname{Var}(z(T)-z(0))
&=\sum_{i=1}^{N}\operatorname{Var}(\varepsilon_i\sqrt{\Delta t})\\
&=N\Delta t=T.
\end{aligned}
$$

독립 정규변수의 합도 정규이므로

$$
\boxed{z(T)-z(0)\sim N(0,\sqrt T),\qquad
\operatorname{SD}(z(T)-z(0))=\sqrt T.}
$$

원문 식 (7)–(9)를 기호 그대로 바로잡으면 다음 세 식이다.

$$
\boxed{\mathbb E[z(T)-z(0)]=0,\qquad
\operatorname{Var}[z(T)-z(0)]=T,\qquad
\operatorname{SD}[z(T)-z(0)]=\sqrt T.}
$$

**독립 증분의 분산은 기간에 비례하고, 표준편차는 기간의 제곱근에 비례한다.** 원문 9쪽에서 긴 기간 $$T$$의 분산을 $$\Delta t$$로 표시한 부분은 $$T$$로 고쳐야 한다.

### 4.3 마르코프 성질의 의미

과거 정보 전체를 $$\mathcal F_t$$라고 하자. 위너 과정에서는 미래 증분이 $$\mathcal F_t$$와 독립이므로 미래 분포를 정할 때 현재 상태만 있으면 된다.

$$
\Pr(z(t+h)\le x\mid\mathcal F_t)
=\Pr(z(t+h)\le x\mid z(t)).
$$

이것이 마르코프 성질이다. 다만 **마르코프 성질만으로 금융시장의 약형 효율성이 증명되는 것은 아니다.** 시장효율성은 이용 가능한 정보, 위험조정 수익률과 거래비용 등을 포함한 경제적 주장이다. 또한 독립 증분은 마르코프 성질보다 강한 조건이다.

## 5. 일반화 위너 과정: 추세와 충격을 더한다

상수 $$a,b$$에 대해

$$
\boxed{dx=a\,dt+b\,dz}
$$

를 생각하자. $$a$$는 단위시간당 기대변화량인 드리프트, $$b^2$$는 단위시간당 분산이다.

원문의 작은 구간 표기로는, 누락된 $$b$$를 복원하여

$$
\boxed{\Delta x=a\Delta t+b\varepsilon\sqrt{\Delta t}},
\qquad \mathbb E[\Delta x]=a\Delta t,
\quad \operatorname{Var}(\Delta x)=b^2\Delta t.
$$

상수계수이므로 구간을 적분하면 정확히

$$
x(t+h)-x(t)=ah+b(z(t+h)-z(t))
=ah+b\sqrt h\,\varepsilon.
$$

따라서

$$
\boxed{x(t+h)-x(t)\sim N(ah,|b|\sqrt h)},
$$

$$
\mathbb E[x(t+h)\mid\mathcal F_t]=x(t)+ah,
\qquad \operatorname{SD}(x(t+h)\mid\mathcal F_t)=|b|\sqrt h.
$$

보통 변동성 계수는 $$b\ge0$$로 놓지만 일반적인 표준편차 식에는 절댓값이 필요하다.

### 수치 예시

$$a=2$$, $$b=3$$, $$h=0.25$$년이면

$$
\Delta x=0.5+1.5\varepsilon,
\qquad \mathbb E[\Delta x]=0.5,
\qquad \operatorname{Var}(\Delta x)=2.25.
$$

드리프트만 존재하면 직선 $$x(t)=x_0+at$$이고, 위너 항을 더하면 그 직선 주변에서 무작위로 움직인다. 장래의 값 자체는 정규분포이므로 음수가 될 가능성도 남는다.

## 6. 주가에는 왜 변화액보다 변화율 모형을 쓰는가?

일반화 위너 과정은 주가가 10이든 100이든 기대변화액 $$a\,dt$$와 충격의 크기 $$b\sqrt{dt}$$가 같다. 하지만 같은 비율의 기대수익률과 변동성을 표현하려면 변화액이 현재 주가에 비례해야 한다.

따라서 상수 $$\mu$$, $$\sigma\ge0$$와 $$S_0>0$$에 대해 다음 모형을 도입한다.

$$
\boxed{dS=\mu S\,dt+\sigma S\,dz}
$$

또는

$$
\boxed{\frac{dS}{S}=\mu\,dt+\sigma\,dz}.
$$

이것이 **기하 브라운 운동(GBM)**이다. 여기서는 배당을 지급하지 않는 주가를 모델링한다.

- \$$\mu$$: 순간적인 기대수익률의 연율 계수.
- \$$\sigma$$: 연율 변동성 계수.
- \$$\sigma^2$$: 순간 수익률의 단위시간당 분산.

### 6.1 변동성이 없을 때

$$\sigma=0$$이면 보통의 미분방정식이 되어

$$
\frac{dS}{S}=\mu\,dt
\quad\Longrightarrow\quad
\int_{S_0}^{S_T}\frac{dS}{S}=\int_0^T\mu\,dt
$$

이므로

$$
\ln\frac{S_T}{S_0}=\mu T,
\qquad S_T=S_0e^{\mu T}.
$$

### 6.2 짧은 구간의 근사

$$\Delta t$$가 작고 구간 초 가격을 $$S$$로 고정하는 Euler 근사는

$$
\frac{\Delta S}{S}\approx\mu\Delta t+\sigma\sqrt{\Delta t}\,\varepsilon_t.
$$

따라서 조건부 근사분포는

$$
\boxed{\frac{\Delta S}{S}\ \Big|\ \mathcal F_t
\ \approx N(\mu\Delta t,\sigma\sqrt{\Delta t}).}
$$

이는 유한기간 단순수익률의 **정확한** 정규분포를 뜻하지 않는다. 너무 큰 시간 간격으로 이 식을 시뮬레이션하면 음의 가격을 만들 수도 있다.

## 7. GBM의 정확한 해: 로그수익률이 정규분포를 따른다 — 보충

확률과정에는 일반 미적분의 연쇄법칙을 그대로 사용할 수 없다. 이토 보조정리를 $$f(S)=\ln S$$에 적용하면

$$
\boxed{d\ln S=\left(\mu-\frac{\sigma^2}{2}\right)dt+\sigma\,dz.}
$$

유도는 [부록 B](#appendix-b)에 있다. 구간을 적분하면

$$
\ln\frac{S_{t+h}}{S}
=\left(\mu-\frac{\sigma^2}{2}\right)h
+\sigma(z(t+h)-z(t)).
$$

따라서

$$
\boxed{S_{t+h}=S\exp\left[
\left(\mu-\frac{\sigma^2}{2}\right)h+\sigma\sqrt h\,\varepsilon\right].}
$$

이는 상수계수 GBM의 격자 시점 사이 **정확한 전이식**이다.

### 7.1 무엇이 정규이고 무엇이 로그정규인가?

여기서 $$\mu$$는 원문 주가모형의 드리프트다. 별도의 $$m$$을 도입하지 않고 로그수익률의 연율 평균을 $$\mu-\sigma^2/2$$로 직접 쓴다. 원문처럼 수익률에는 $$r$$를 쓰되, 아래 비교표와 비교식에서만 설명용 윗첨자 $$\mathrm{log}$$와 $$\mathrm{simple}$$로 종류를 구별한다.

| 대상 | 정확한 성질 |
|---|---|
| 로그수익률 $$r^{\mathrm{log}}_{t,t+h}=\ln(S_{t+h}/S)$$ | $$N((\mu-\sigma^2/2)h,\sigma\sqrt h)$$ |
| 총수익배율 $$S_{t+h}/S$$ | 로그정규분포 |
| 주가 $$S_{t+h}$$, 현재 가격이 주어졌을 때 | 로그정규분포, 항상 양수 |
| 단순수익률 $$r^{\mathrm{simple}}_{t,t+h}=S_{t+h}/S-1$$ | 로그정규변수에서 1을 뺀 분포, 하한은 -1 |

**주가가 정규분포라는 말과 로그수익률이 정규분포라는 말은 다르다.**

![일반화 위너 과정과 기하 브라운 운동]({{ '/assets/img/investments/processes.png' | relative_url }})
_같은 정규 충격으로 생성한 모의 경로다. 왼쪽은 가산적 과정, 오른쪽은 승산적 GBM이다. GBM은 정확한 전이식을 사용했다. 두 패널의 세로축 단위는 서로 다르다._

### 7.2 기대수익률과 기대 로그수익률은 다르다

정규변수 $$X\sim N(a,\sqrt v)$$에 대해 $$\mathbb E[e^X]=e^{a+v/2}$$이므로

$$
\mathbb E[S_{t+h}\mid\mathcal F_t]=Se^{\mu h},
$$

$$
\mathbb E[r^{\mathrm{simple}}_{t,t+h}]=e^{\mu h}-1,
\qquad \mathbb E[r^{\mathrm{log}}_{t,t+h}]=\left(\mu-\frac{\sigma^2}{2}\right)h.
$$

즉, $$\mu=12\%$$는 1년 단순수익률의 기대값이 정확히 12%라는 의미가 아니다. 이 모형에서 그 값은 $$e^{0.12}-1\approx12.75\%$$다.

또한

$$
\operatorname{Var}(S_{t+h}\mid\mathcal F_t)
=S^2e^{2\mu h}(e^{\sigma^2h}-1).
$$

평균·분산의 상세 계산은 [부록 C](#appendix-c)에 있다.

## 8. 연율을 일별 수치로 바꾸기

한 해를 252거래일로 놓고, 독립·동일분포인 로그증분을 가정하자. $$\mu$$를 앞 절의 **주가모형 드리프트**로 유지하면 일별 로그수익률 $$r_{\mathrm{day}}$$은

$$
r_{\mathrm{day}}=\ln\frac{S_{t+1/252}}S
\sim N\left(\frac{\mu-\sigma^2/2}{252},\frac{\sigma}{\sqrt{252}}\right).
$$

따라서

$$
\boxed{\mathbb E[r_{\mathrm{day}}]=\frac{\mu-\sigma^2/2}{252},\qquad
\operatorname{SD}(r_{\mathrm{day}})=\frac\sigma{\sqrt{252}}.}
$$

원문 13쪽은 연율 **로그수익률 평균**에도 $$\mu$$라는 기호를 사용하므로 그 문단의 정의만 따르면 일별 평균을 $$\mu/252$$로 쓴다. 그러나 원문 11–12쪽 주가모형의 $$\mu$$를 이어서 사용할 때에는 위처럼 $$\sigma^2/2$$ 보정을 해야 한다. **같은 문자라도 두 문단의 정의가 같지 않다.**

일별 로그수익률의 평균을 $$\mathbb E[r_i]$$, 분산을 $$\operatorname{Var}(r_i)$$로 직접 표시하면, $$n$$개 독립·동일분포 증분에 대해

$$
\mathbb E\left[\sum_{i=1}^nr_i\right]=n\mathbb E[r_i],
\qquad
\operatorname{Var}\left(\sum_{i=1}^nr_i\right)=n\operatorname{Var}(r_i).
$$

평균은 $$n$$배, 표준편차는 $$\sqrt n$$배다. 평균의 선형 연율화와 실현 복리수익률은 다르다. **단순수익률**을 $$r_i$$로 정의한 경우 실현 복리수익률은 $$\prod_i(1+r_i)-1$$로 계산한다. 자기상관이 있으면 분산의 제곱근 환산도 수정해야 한다.

## 9. 정규분포가 극단적 움직임을 과소평가하는 이유

### 9.1 강의의 급락 계산을 일관되게 다시 계산하기

원문은 연율 12%, 변동성 15%를 사용하는 정규근사 예시를 제시한다. 여기서는 **강의의 계산방식을 재현하는 가상 일별 단순수익률 모형**을 먼저 정의한다.

$$
r\sim N\left(\frac{0.12}{252},\frac{0.15}{\sqrt{252}}\right).
$$

이 수치는 현재 시장 추정치가 아니다. 평균과 표준편차는

$$
\mu_d\approx0.00047619,
\qquad \sigma_d\approx0.00944911.
$$

원문처럼 $$X=(r-0.00047619)/0.00944911$$로 표준화한다. 그러면 $$X\sim N(0,1)$$이고, 23% 이상 하락할 확률은

$$
\begin{aligned}
\Pr(r\le-0.23)
&=N\left(\frac{-0.23-0.00047619}{0.00944911}\right)\\
&=N(-24.3913)\\
&\approx1.06\times10^{-131}.
\end{aligned}
$$

이 값은 $$10^{-127}$$이 아니다. 원문의 부호·환산식·표준화 숫자가 서로 일치하지 않아 다시 계산했다. 또한 원문의 역사적 사건에 대한 지수명과 하락률의 조합은 여기서 확정하지 않고 **-23%라는 가상 임계값**으로만 사용한다.

같은 가상 정규모형에서 14% 이상 상승할 확률은

$$
\Pr(r\ge0.14)=1-N(14.7658)
\approx1.22\times10^{-49}.
$$

원문의 두 역사적 사례는 서로 다른 지수에 관한 것이므로, 이것을 Nasdaq 고유의 모수를 추정한 결과로 해석하면 안 된다.

### 9.2 GBM을 정확하게 적용하면 임계값도 로그로 바꿔야 한다

이번에는 $$\mu=0.12,\sigma=0.15$$인 GBM을 가정하자. -23%의 **단순수익률**은

$$
r^{\mathrm{log}}\le\ln(0.77)\approx-0.261365
$$

와 같은 사건이다. 따라서

$$
\Pr(r\le-0.23)
=N\left(
\frac{\ln(0.77)-(0.12-0.15^2/2)/252}{0.15/\sqrt{252}}
\right)
\approx2.96\times10^{-169}.
$$

두 확률의 차이는 계산오류가 아니라 **서로 다른 분포 가정** 때문이다. 어느 계산에서도 고정 변동성 정규모형은 극단사건에 매우 작은 확률을 부여한다.

### 9.3 두꺼운 꼬리

정규분포도 극단적 사건에 양의 확률을 준다. 문제는 실제 금융수익률의 큰 변동을 설명하기에 그 확률이 지나치게 작을 수 있다는 것이다.

![정규분포와 두꺼운 꼬리 분포]({{ '/assets/img/investments/fat-tails.png' | relative_url }})
_평균 0·분산 1로 맞춘 정규분포와 자유도 5의 Student t 분포를 비교했다. 오른쪽은 세로축을 로그척도로 표시하여 꼬리 차이를 드러낸다. 실제 수익률에 적합한 결과가 아닌 이론 비교다._

정규모형의 한계는 점프, 시간에 따라 달라지는 변동성, 비대칭성 등과 연결된다. 원문의 핵심은 정규분포를 무조건 폐기하는 것이 아니라, **모형이 주는 작은 확률을 현실의 불가능성과 동일시하지 말라**는 것이다.

## 10. 실제 데이터를 분석하기 전: 무엇을 관측했는가?

가격수준 $$P_t$$와 수익률 $$r_t$$는 서로 다른 시계열이다. 보통 수익률을 분석할 때에는

$$
r_t=\frac{P_t-P_{t-1}+D_t}{P_{t-1}},
\qquad r_t^{\mathrm{log}}=\ln(1+r_t)
$$

처럼 변환한다. 이 비교식에서 $$r_t$$는 단순수익률이며, 로그수익률에만 설명용 윗첨자를 붙였다. 배당·분할을 일관되게 조정한 총수익지수를 이용할 수도 있다.

원문은 추세를 확인하고 제거하라고 설명한다. 이를 “모든 가격자료에 선형추세를 회귀해서 빼면 된다”는 뜻으로 읽어서는 안 된다. 로그차분, 결정적 추세 제거, 계절조정 등은 자료의 생성과정에 따라 구별된다.

### 정상성, 에르고드성, i.i.d.

| 가정 | 의미 | 분석에서의 역할 |
|---|---|---|
| 약정상성 | 유한한 평균·분산이 시간에 따라 일정하고 공분산은 시차에만 의존 | 시차별 의존구조를 안정적으로 요약 |
| 강정상성 | 유한차원 결합분포가 시간 이동에 불변 | 전체 확률법칙의 시간 불변성 |
| 에르고드성 | 적절한 조건에서 시간평균으로 모집단 특성을 회복 가능 | 하나의 긴 시계열로 평균 등을 추정 |
| i.i.d. | 관측치가 서로 독립이고 동일한 분포를 가짐 | 단순한 평균·분산·표준오차 공식 |

정상성만으로 독립성이 보장되지는 않는다. 구조변화가 있다면 관측기간을 늘리는 것이 서로 다른 분포를 섞는 결과가 될 수도 있다.

## 11. 경험분포와 히스토그램

관측 수익률을 $$r_1,\ldots,r_N$$이라 하자. 여기서 $$r_i$$는 단순수익률 또는 로그수익률 중 **하나를 정해 일관되게 사용한 자료**다.

경험누적분포함수는

$$
\boxed{\widehat F_N(x)=\frac1N\sum_{i=1}^{N}\mathbf1\{r_i\le x\}.}
$$

히스토그램을 위해 $$\underline{x}<\overline{x}$$인 관측범위를 $$K$$등분하면 구간 폭은

$$
\Delta x=\frac{\overline{x}-\underline{x}}K.
$$

경계값이 두 구간에 중복 집계되지 않도록 마지막을 제외한 구간은 왼쪽을 포함하고 오른쪽을 제외한다. $$k$$번째 구간의 관측 수가 $$N_k$$이면

$$
\widehat p_k=\frac{N_k}{N},\qquad \sum_k\widehat p_k=1.
$$

확률밀도 곡선과 비교하려면 막대 높이를

$$
\widehat f_k=\frac{N_k}{N\Delta x}
$$

로 정해야 한다. 그러면

$$
\sum_k\widehat f_k\Delta x=1.
$$

**구간 확률과 확률밀도는 단위가 다르다.** $$N_k/N$$를 그대로 정규밀도와 비교하면 잘못된 그래프가 될 수 있다. 모든 관측값이 같다면 위 등분 폭은 0이므로 별도의 구간을 정해야 한다.

## 12. 표본평균·분산·왜도·첨도

### 12.1 평균

원문 식 (33)처럼 **표본평균**을 $$\mu$$로 쓴다.

$$
\boxed{\mu=\frac1N\sum_{i=1}^{N}r_i.}
$$

이 절 이후의 표본통계 설명에서 $$\mu$$는 주가모형의 드리프트가 아니라 위 식으로 계산한 표본평균이다. 표본평균은 자료가 바뀌면 값도 바뀌는 확률변수다. 증명에서 모집단 평균과 비교할 때만 보충 기호 $$\mu_{\mathrm{pop}}$$을 사용한다.

### 12.2 분산: 원문의 N과 불편보정의 N-1

원문 식 (34)의 경험분산 표기를 그대로 유지한다.

$$
\boxed{\sigma^2=\frac1N\sum_{i=1}^{N}(r_i-\mu)^2.}
$$

여기서 $$\sigma$$는 위 경험분산의 제곱근이다. 모집단 분산과 구별해야 하는 증명에서는 모집단 분산을 보충 기호 $$\sigma_{\mathrm{pop}}^2$$로 쓴다.

i.i.d. 표본에서 모집단 분산의 불편추정량은 $$N>1$$일 때

$$
s^2=\frac1{N-1}\sum_{i=1}^{N}(r_i-\mu)^2
=\frac{N}{N-1}\sigma^2.
$$

$$s^2$$는 **원문의 식을 바꾸기 위한 기호가 아니라**, 불편분산을 추가 설명하기 위한 기호다. 원문의 $$1/N$$은 경험분산으로 올바르므로 그대로 두었다. $$N-1$$이 필요한 이유는 [부록 D](#appendix-d)에서 증명한다.

### 12.3 왜도와 첨도

원문의 이름 $$\mathrm{skew}$$와 $$\mathrm{kurt}$$를 그대로 쓰고, 잘못된 분모만 수정한다. $$\sigma>0$$일 때

$$
\boxed{\mathrm{skew}=\frac{\frac1N\sum_{i=1}^{N}(r_i-\mu)^3}{\sigma^3}},
$$

$$
\boxed{\mathrm{kurt}=\frac{\frac1N\sum_{i=1}^{N}(r_i-\mu)^4}{\sigma^4}}.
$$

왜도의 분모는 $$\sqrt{\sigma^3}$$가 아니라 $$\sigma^3$$이고, 첨도의 분모는 $$\sigma^2$$가 아니라 $$\sigma^4$$다. 분자와 분모의 단위가 같아져 두 통계량이 무차원이 되어야 한다는 점으로도 확인할 수 있다.

초과첨도는 $$\mathrm{kurt}-3$$이다.

| 지표 | 해석 |
|---|---|
| 음의 왜도 | 왼쪽 극단값이 3차 중심적률에 더 크게 기여 |
| 양의 왜도 | 오른쪽 극단값이 3차 중심적률에 더 크게 기여 |
| 정규분포의 모집단 첨도 | 3 |
| 정규분포의 모집단 초과첨도 | 0 |
| 큰 첨도 | 평균에서 멀리 떨어진 관측치의 4차 적률 기여가 큼 |

첨도를 단지 “봉우리의 높이”로 해석하면 부정확하다. 왜도 0이 반드시 대칭성을 뜻하지 않고, 첨도 3이 정규분포를 보장하지도 않는다. 위 식은 적률형 표본통계량으로, 유한표본 불편보정을 적용한 왜도·첨도와는 다르다.

## 13. 표준오차: 추정한 평균은 얼마나 정확한가?

$$r_i$$가 i.i.d.이고 $$\mathbb E[r_i]=\mu_{\mathrm{pop}}$$, $$\operatorname{Var}(r_i)=\sigma_{\mathrm{pop}}^2<\infty$$라고 하자.

기대값은

$$
\mathbb E[\mu]=\frac1N\sum_i\mathbb E[r_i]=\mu_{\mathrm{pop}}.
$$

분산은 독립성을 이용하여

$$
\begin{aligned}
\operatorname{Var}(\mu)
&=\operatorname{Var}\left(\frac1N\sum_ir_i\right)\\
&=\frac1{N^2}\sum_i\sigma_{\mathrm{pop}}^2
=\frac{\sigma_{\mathrm{pop}}^2}{N}.
\end{aligned}
$$

따라서

$$
\boxed{\operatorname{SE}(\mu)=\frac{\sigma_{\mathrm{pop}}}{\sqrt N},
\qquad \widehat{\operatorname{SE}}(\mu)=\frac{s}{\sqrt N}.}
$$

원문의 경험표준편차 $$\sigma$$를 대입한 큰 표본 근사는 $$\sigma/\sqrt N$$이다. 불편분산으로 보정하면 위 식의 $$s/\sqrt N$$을 사용하며, $$s=\sqrt{N/(N-1)}\,\sigma$$다.

**표준편차는 수익률 자체의 흩어짐이고, 표준오차는 평균 추정치의 흩어짐이다.** 표본 수가 4배가 되면 평균의 표준오차는 절반이 되지만 개별 수익률의 변동성이 절반이 되는 것은 아니다.

![표본크기와 평균의 표준오차]({{ '/assets/img/investments/standard-error.png' | relative_url }})
_일별 표준편차를 1%로 고정한 i.i.d. 모형. 관측 수를 늘릴수록 표준오차는 감소하지만 감소 속도는 점차 느려진다._

### 수치 예시와 신뢰구간

$$N=252$$, $$\mu=0.0005$$, $$s=0.01$$이면

$$
\widehat{\operatorname{SE}}(\mu)=\frac{0.01}{\sqrt{252}}
\approx0.00062994.
$$

이는 0.062994%포인트다. 큰 표본에서 정규근사한 95% 신뢰구간은

$$
0.0005\pm1.96(0.00062994)
\approx[-0.0007347,\;0.0017347].
$$

백분율로는 약 $$[-0.0735\%,0.1735\%]$$다. 정규 i.i.d. 표본에서 분산이 미지이면 정확한 구간은 $$1.96$$ 대신 $$t_{N-1,0.975}$$를 사용한다. 이 신뢰구간은 **모집단 평균**을 추정하는 구간이지 내일 수익률의 예측구간이 아니다.

### 더 자주 관측하면 연율 평균도 정확해지는가? — 보충

고정된 기간 $$T$$년의 GBM을 $$N$$등분해 로그수익률 $$r_i$$를 관측하자. 이 보충의 $$\mu,\sigma$$는 앞의 표본통계량이 아니라 GBM의 모수로 다시 사용한다. $$\widehat\mu_{\mathrm{log}}$$는 연율 로그평균 $$\mu-\sigma^2/2$$의 추정량이다. 그러면

$$
\widehat\mu_{\mathrm{log}}=\frac1T\sum_{i=1}^{N}r_i
=\frac{\ln(S_T/S_0)}T,
\qquad \operatorname{SE}(\widehat\mu_{\mathrm{log}})=\frac\sigma{\sqrt T}.
$$

고정된 $$T$$ 안에서 관측 빈도만 높여도 연율 드리프트의 표준오차가 줄어드는 것은 아니다. **같은 간격의 관측을 더 오래 모으는 것**과 **같은 기간을 더 잘게 나누는 것**을 구별해야 한다. [부록 E](#appendix-e)에 전개했다.

## 14. 조건부 분포: 평균과 변동성이 시간에 따라 바뀐다

전체 기간에 하나의 평균과 분산을 적용하는 분석은 무조건부 분포를 요약한다. 과거 정보 $$\mathcal F_{t-1}$$를 이용해

$$
r_t=\mu_t+\sigma_t\varepsilon_t,
\qquad
\mathbb E[\varepsilon_t\mid\mathcal F_{t-1}]=0,
\quad \operatorname{Var}(\varepsilon_t\mid\mathcal F_{t-1})=1
$$

로 쓰면

$$
\mathbb E[r_t\mid\mathcal F_{t-1}]=\mu_t,
\qquad \operatorname{Var}(r_t\mid\mathcal F_{t-1})=\sigma_t^2.
$$

원문은 월별 표본평균·표준편차를 계산하는 간단한 접근을 소개한다. 다만 같은 달의 모든 수익률로 계산한 평균·변동성은 **사후 요약치**다. 그 달이 시작되기 전 알 수 있었던 예측치와는 다르다.

과거 $$M$$개 관측만 사용하는 예측용 이동창 예시는

$$
\widehat\mu_t=\frac1M\sum_{j=1}^{M}r_{t-j},
$$

$$
\widehat\sigma_t^2=\frac1{M-1}\sum_{j=1}^{M}(r_{t-j}-\widehat\mu_t)^2.
$$

미래 데이터를 섞지 않는 대신, 창이 짧으면 추정이 불안정하고 길면 변화에 늦게 반응한다.

### 14.1 변동성 군집과 비대칭성

원문은 높은 변동성이 이어지는 현상과 하락기에 변동성이 높아지는 경향을 소개한다. 수익률의 자기상관이 작아도 $$r_t^2$$나 $$|r_t|$$에는 시계열 의존성이 남을 수 있다. 이런 의존성은 단순 i.i.d. 가정을 점검할 이유가 된다.

### 14.2 조건부 정규인데 무조건부로는 꼬리가 두꺼울 수 있다 — 보충

$$r=\sqrt V\,\varepsilon$$이고 $$V>0$$, $$\varepsilon\sim N(0,1)$$가 독립이라고 하자. 각 $$V$$가 주어졌을 때는 정규분포지만 $$V$$ 자체가 변하면 여러 분산의 정규분포가 섞인다.

$$
\mathbb E[r^2]=\mathbb E[V],\qquad
\mathbb E[r^4]=3\mathbb E[V^2].
$$

따라서 유한한 관련 적률이 존재할 때 무조건부 첨도는

$$
\frac{\mathbb E[r^4]}{(\mathbb E[r^2])^2}
=3\frac{\mathbb E[V^2]}{(\mathbb E[V])^2}
=3\left(1+\frac{\operatorname{Var}(V)}{(\mathbb E[V])^2}\right)\ge3.
$$

분산이 실제로 변하면 부등호가 엄격하다. 이는 시변 변동성이 정규분포보다 큰 첨도를 만들어낼 수 있음을 보여준다. 모든 두꺼운 꼬리 현상을 이 한 가지 모형으로 설명한다는 뜻은 아니다.

## 15. 자기상관이 있으면 표준오차도 달라진다 — 보충

약정상 수익률의 시차 $$k$$ 자기공분산을 $$\gamma_k=\operatorname{Cov}(r_t,r_{t-k})$$라 하자. 그러면

$$
\boxed{\operatorname{Var}(\mu)
=\frac1{N^2}\left[N\gamma_0+2\sum_{k=1}^{N-1}(N-k)\gamma_k\right].}
$$

i.i.d.이면 $$k\ge1$$에서 $$\gamma_k=0$$이어서 익숙한 $$\gamma_0/N$$이 된다. 자기상관이 있으면 단순한 $$s/\sqrt N$$이 정확하지 않을 수 있다. 합을 전개한 증명은 [부록 F](#appendix-f)에 있다.

## 16. 전체 연결

| 단계 | 얻는 것 | 확인할 가정 |
|---|---|---|
| 수익률 정의 | 비교 가능한 투자성과 | 배당·분할·복리 처리 |
| 확률분포 | 미래 결과의 가능성 | 분포와 꼬리 형태 |
| 위너 과정 | 시간에 따른 독립 정규증분 | 분산의 시간 비례 |
| GBM | 양의 주가와 정규 로그수익률 | 상수 드리프트·변동성 |
| 경험분포·표본적률 | 자료에 나타난 분포 요약 | 표본범위·정상성 |
| 표준오차 | 추정치의 불확실성 | 독립성·표본길이 |
| 조건부 분석 | 시점별 분포 변화 | 가용 정보와 예측 시점 |

투자에서는 수익률의 크기뿐 아니라 **분포 가정이 타당한지, 모수 추정이 얼마나 불확실한지**도 함께 판단해야 한다.

## 원문 정오표와 해석상 주의

쪽수는 첨부된 27쪽 PDF 기준이다.

| 위치 | 원문 표기·서술 | 이 글의 수정·구별 |
|---|---|---|
| 4쪽 | 정규분포 2표준편차 확률 95.54% | 약 95.45% |
| 9쪽, 식 (8)·(9) | 긴 기간 $$T$$ 증분의 표준편차·분산에 $$\sqrt{\Delta t},\Delta t$$ | $$\sqrt T,T$$ |
| 9쪽, 식 (11) 앞 | 확률항을 포함한다고 하면서 결정적 식 제시 | 확률항을 제거한 경우로 해석 |
| 9쪽, 식 (14) | 충격항에서 $$b$$ 누락 | $$\Delta x=a\Delta t+b\varepsilon\sqrt{\Delta t}$$ |
| 11쪽, 식 (26) | $$dSS$$ 형태의 오기 | $$dS/S$$ |
| 12쪽, 식 (29) | 정규분포의 두 번째 인수에 표준편차 사용 | 오류가 아닌 표기 관례이므로 유지. 이 글은 $$N(\text{평균},\text{표준편차})$$ 사용 |
| 13쪽 | 로그수익률 평균과 가격과정 드리프트에 같은 기호 사용 | GBM의 $$\mu$$를 유지할 때 로그평균은 $$\mu-\sigma^2/2$$라고 직접 표시 |
| 15쪽 | 일별 표준편차를 $$0.15\sqrt{252}$$로 기재 | $$0.15/\sqrt{252}$$ |
| 15쪽 | 급락 사건을 $$r<0.23$$로 기재 | $$r<-0.23$$ |
| 15쪽, 식 (31)·(32) | 표준화 값 -23 및 $$10^{-127}$$ | 제시된 정규근사 모수로 -24.3913, 약 $$1.06\times10^{-131}$$ |
| 17쪽 | 수익률 그래프 세로축에는 weekly, 캡션에는 daily | 원자료 없이 빈도를 단정하지 않고 재사용하지 않음 |
| 20쪽 | 히스토그램 구간 폭에 최댓값·최솟값의 합 | $$(\overline{x}-\underline{x})/K$$, 경계 중복 방지 |
| 22쪽, 식 (34) | 분산의 분모 $$N$$ | 경험분산으로 유지하고 불편분산 $$N-1$$과 구분 |
| 22쪽, 식 (35) | 왜도 분모 $$\sqrt{\sigma^3}$$ | 표준편차 기준 $$\sigma^3$$ |
| 22쪽, 식 (36) | 첨도 분모 $$\sigma^2$$ | 표준편차 기준 $$\sigma^4$$ |

## 자료와 보충 범위

- 기본 자료: 업로드된 「투자론 2장.pdf」, MIT 15.433 Investments, Class 2, Spring 2003. 확률분포 3–5쪽, 확률과정 6–13쪽, 극단사건 14–16쪽, 경험분포와 통계적 추정 17–25쪽을 중심으로 정리했다.
- 26–27쪽의 요약과 다음 장 연결은 반영하고, 강의자료에 본문이 없는 별도 교재 문항·페이지 목록은 제외했다.
- GBM의 정규 로그수익률과 로그정규 주가 관계는 [Columbia University의 GBM 강의노트](https://www.columbia.edu/~ks20/FE-Notes/4700-07-Notes-GBM.pdf)와 교차 확인했다. 그 자료의 로그 드리프트는 이 글의 GBM 계수로 쓰면 \$$\mu-\sigma^2/2$$에 해당하므로 정의를 구분했다.
- 모든 수치 예시는 설명용 가정이며 현재 시장의 추정치가 아니다. 그래프는 직접 계산한 이론 그래프 또는 모의실험이다.

## 이 장의 수학 부록
{: #chapter-appendix }

본문에서 분리한 긴 증명과 보충 유도다.

### 부록 B. 이토 보조정리로 GBM의 해 구하기
{: #appendix-b }

#### B.1 출발점과 목표

상수계수 주가모형

$$
dS=\mu S\,dt+\sigma S\,dz,\qquad S_0>0
$$

에서 $$S$$의 명시적 표현을 얻고 싶다. $$f(S)=\ln S$$를 적용하면 곱셈적 변화를 덧셈적 변화로 바꿀 수 있다.

#### B.2 일반 연쇄법칙과 다른 항

이토 보조정리는 충분히 매끄러운 $$f(t,S)$$에 대해

$$
df=\left(f_t+\mu Sf_S+\frac12\sigma^2S^2f_{SS}\right)dt
+\sigma Sf_S\,dz
$$

라고 말한다. 이를 정리로 사용하여 GBM의 해를 유도한다. 이토 보조정리 자체의 확률적분론적 증명은 이 부록의 범위 밖이다.

왜 2차 미분항이 남는지는 테일러 전개로 직관을 얻을 수 있다.

$$
\Delta f\approx f_S\Delta S+\frac12f_{SS}(\Delta S)^2.
$$

위너 증분은 크기가 $$\sqrt{\Delta t}$$ 수준이므로 그 제곱은 $$\Delta t$$ 수준이다. 극한에서 이 2차항을 무시할 수 없다. 흔히 쓰는

$$
(dz)^2=dt,\qquad dt\,dz=0,\qquad(dt)^2=0
$$

은 보통의 실수 미분을 곱하는 등식이 아니라 이토 계산과 이차변동을 요약하는 표기다.

#### B.3 로그에 적용

$$f(S)=\ln S$$이면

$$
f_S=\frac1S,\qquad f_{SS}=-\frac1{S^2},\qquad f_t=0.
$$

대입하면

$$
\begin{aligned}
d\ln S
&=\left[\frac1{S}\mu S
+\frac12\left(-\frac1{S^2}\right)\sigma^2S^2\right]dt
+\frac1{S}\sigma S\,dz\\
&=\left(\mu-\frac12\sigma^2\right)dt+\sigma\,dz.
\end{aligned}
$$

0부터 $$T$$까지 적분하면

$$
\ln S_T-\ln S_0
=\left(\mu-\frac12\sigma^2\right)T+\sigma z(T).
$$

따라서

$$
\boxed{S_T=S_0\exp\left[\left(\mu-\frac12\sigma^2\right)T+\sigma z(T)\right].}
$$

지수함수는 양수이므로 $$S_T>0$$이다. 엄밀하게는 위 양의 후보해에 이토 보조정리를 적용해 원래 확률미분방정식을 만족함을 확인할 수도 있다.

#### B.4 로그 드리프트의 보정항

$$-\sigma^2/2$$는 임의로 추가한 손실항이 아니다. 로그함수의 오목성과 확률적 변동이 결합해 생기는 이토 보정이다.

$$
\mathbb E\left[\ln\frac{S_T}{S_0}\right]
=\left(\mu-\frac12\sigma^2\right)T,
$$

반면

$$
\ln\mathbb E\left[\frac{S_T}{S_0}\right]=\mu T.
$$

둘의 차이는 $$\sigma^2T/2$$이다. 로그의 기대값과 기대값의 로그는 같지 않다.

### 부록 C. 정규변수의 지수적률과 GBM의 평균·분산
{: #appendix-c }

#### C.1 정규변수의 적률생성함수

이 부록에서 $$X$$는 표준화 변수가 아닌 보조 정규변수다. $$X\sim N(a,\sqrt v)$$이고 분산 $$v>0$$이면

$$
\mathbb E[e^{uX}]
=\frac1{\sqrt{2\pi v}}\int_{-\infty}^{\infty}
\exp\left[uy-\frac{(y-a)^2}{2v}\right]dy.
$$

지수 부분을 완전제곱으로 정리한다.

$$
uy-\frac{(y-a)^2}{2v}
=-\frac{(y-a-uv)^2}{2v}+ua+\frac12u^2v.
$$

따라서

$$
\mathbb E[e^{uX}]
=e^{ua+u^2v/2}\frac1{\sqrt{2\pi v}}
\int_{-\infty}^{\infty}\exp\left[-\frac{(y-a-uv)^2}{2v}\right]dy.
$$

뒤의 정규밀도 적분은 1이므로

$$
\boxed{\mathbb E[e^{uX}]=e^{ua+u^2v/2}.}
$$

$$v=0$$인 경우에도 퇴화확률변수를 직접 계산하면 같은 식이 성립한다.

#### C.2 주가의 1차 적률

$$X=\ln(S_T/S_0)$$라 하면

$$
X\sim N\left(\left(\mu-\frac12\sigma^2\right)T,\sigma\sqrt T\right).
$$

$$u=1$$을 대입하여

$$
\begin{aligned}
\mathbb E[S_T]
&=S_0\mathbb E[e^X]\\
&=S_0\exp\left[\left(\mu-\frac12\sigma^2\right)T+\frac12\sigma^2T\right]\\
&=S_0e^{\mu T}.
\end{aligned}
$$

#### C.3 주가의 2차 적률과 분산

$$u=2$$를 대입하면

$$
\begin{aligned}
\mathbb E[S_T^2]
&=S_0^2\mathbb E[e^{2X}]\\
&=S_0^2\exp\left[2\left(\mu-\frac12\sigma^2\right)T+2\sigma^2T\right]\\
&=S_0^2e^{(2\mu+\sigma^2)T}.
\end{aligned}
$$

따라서

$$
\begin{aligned}
\operatorname{Var}(S_T)
&=\mathbb E[S_T^2]-\mathbb E[S_T]^2\\
&=S_0^2e^{(2\mu+\sigma^2)T}-S_0^2e^{2\mu T}\\
&=\boxed{S_0^2e^{2\mu T}(e^{\sigma^2T}-1)}.
\end{aligned}
$$

단순수익률 $$r=S_T/S_0-1$$은

$$
\mathbb E[r]=e^{\mu T}-1,
\qquad \operatorname{Var}(r)=e^{2\mu T}(e^{\sigma^2T}-1).
$$

작은 $$T$$에 대해서는 $$e^{\mu T}-1\approx\mu T$$, $$e^{\sigma^2T}-1\approx\sigma^2T$$이므로 2장의 단기 정규근사에서 사용한 평균·분산과 1차 수준에서 연결된다.

### 부록 D. 표본분산의 분모가 N-1인 이유
{: #appendix-d }

$$r_1,\ldots,r_N$$이 i.i.d.이고 평균 $$\mu_{\mathrm{pop}}$$, 유한한 분산 $$\sigma_{\mathrm{pop}}^2$$를 갖는다고 하자.

우선 항등식

$$
\boxed{\sum_{i=1}^{N}(r_i-\mu)^2
=\sum_{i=1}^{N}(r_i-\mu_{\mathrm{pop}})^2-N(\mu-\mu_{\mathrm{pop}})^2}
$$

을 보이자. $$r_i-\mu=(r_i-\mu_{\mathrm{pop}})-(\mu-\mu_{\mathrm{pop}})$$이므로

$$
\begin{aligned}
\sum_i(r_i-\mu)^2
&=\sum_i(r_i-\mu_{\mathrm{pop}})^2
-2(\mu-\mu_{\mathrm{pop}})\sum_i(r_i-\mu_{\mathrm{pop}})
+N(\mu-\mu_{\mathrm{pop}})^2\\
&=\sum_i(r_i-\mu_{\mathrm{pop}})^2
-2N(\mu-\mu_{\mathrm{pop}})^2+N(\mu-\mu_{\mathrm{pop}})^2.
\end{aligned}
$$

마지막 줄에서 $$\sum_i(r_i-\mu_{\mathrm{pop}})=N(\mu-\mu_{\mathrm{pop}})$$를 사용했다.

양변의 기대값을 취하면

$$
\begin{aligned}
\mathbb E\left[\sum_i(r_i-\mu)^2\right]
&=N\sigma_{\mathrm{pop}}^2-N\mathbb E[(\mu-\mu_{\mathrm{pop}})^2]\\
&=N\sigma_{\mathrm{pop}}^2-N\operatorname{Var}(\mu)\\
&=N\sigma_{\mathrm{pop}}^2-\sigma_{\mathrm{pop}}^2\\
&=(N-1)\sigma_{\mathrm{pop}}^2.
\end{aligned}
$$

따라서

$$
\mathbb E\left[\frac1{N-1}\sum_i(r_i-\mu)^2\right]=\sigma_{\mathrm{pop}}^2.
$$

분모를 $$N$$으로 나누면 기대값은 $$(N-1)\sigma_{\mathrm{pop}}^2/N$$이 되어 모집단 분산을 평균적으로 작게 추정한다. 반대로 모집단 평균을 이미 알아서 $$\mu$$ 대신 $$\mu_{\mathrm{pop}}$$를 쓴다면 분모 $$N$$이 맞는다.

이 불편성에는 정규성은 필요하지 않다. 다만 $$s^2$$가 불편이라고 해서 $$s=\sqrt{s^2}$$도 자동으로 불편인 것은 아니다.

### 부록 E. 평균의 표준오차와 관측기간
{: #appendix-e }

#### E.1 고정된 간격의 관측을 늘리는 경우

독립이며 분산이 $$\sigma_{\mathrm{pop}}^2$$로 같은 수익률 $$N$$개에 대해

$$
\operatorname{Var}\left(\frac1N\sum_ir_i\right)
=\frac1{N^2}\sum_i\sigma_{\mathrm{pop}}^2=\frac{\sigma_{\mathrm{pop}}^2}{N}.
$$

$$N$$을 4배로 늘리면 표준오차는 절반이다. 이때 관측 간격은 같고 총 관측기간이 4배로 늘어난다.

#### E.2 고정된 기간을 더 잘게 관측하는 경우

$$d\ln S=\left(\mu-\frac12\sigma^2\right)dt+\sigma\,dz$$에서 총 관측기간을 $$T$$, 관측 간격을 $$h=T/N$$으로 두자. 각 로그수익률은

$$
r_i=\left(\mu-\frac12\sigma^2\right)h+\sigma\sqrt h\,\varepsilon_i.
$$

이 부록에서 $$\mu$$는 다시 GBM의 드리프트다. 로그수익률 $$r_i$$의 연율 평균 $$\mu-\sigma^2/2$$에 대한 추정량을 보충 기호 $$\widehat\mu_{\mathrm{log}}$$로 쓰면

$$
\widehat\mu_{\mathrm{log}}=\frac{\frac1N\sum_i r_i}{h}
=\frac1{Nh}\sum_ir_i
=\frac1T\sum_ir_i.
$$

그 분산은

$$
\begin{aligned}
\operatorname{Var}(\widehat\mu_{\mathrm{log}})
&=\frac1{T^2}\sum_i\operatorname{Var}(r_i)\\
&=\frac1{T^2}N\sigma^2h\\
&=\frac{\sigma^2}{T}.
\end{aligned}
$$

따라서

$$
\boxed{\operatorname{SE}(\widehat\mu_{\mathrm{log}})=\frac\sigma{\sqrt T}.}
$$

$$T$$가 고정되어 있으면 $$N$$을 늘려도 이 표준오차는 변하지 않는다. 변동성 추정은 다른 문제이며, 실제 고빈도 자료에서는 시장미시구조 잡음도 고려해야 한다.

### 부록 F. 자기상관이 있을 때 표본평균의 분산
{: #appendix-f }

약정상 시계열에서 $$\gamma_k=\operatorname{Cov}(r_t,r_{t-k})$$라 하자. 평균의 분산을 이중합으로 쓰면

$$
\operatorname{Var}(\mu)
=\frac1{N^2}\sum_{i=1}^{N}\sum_{j=1}^{N}\operatorname{Cov}(r_i,r_j).
$$

대각선 $$i=j$$에는 $$\gamma_0$$가 $$N$$번 나타난다. 시차가 $$k$$인 쌍 $$j=i+k$$는 $$N-k$$개이며 반대 방향까지 포함하면 $$2(N-k)$$개다.

그러므로

$$
\boxed{\operatorname{Var}(\mu)
=\frac1{N^2}\left[N\gamma_0+2\sum_{k=1}^{N-1}(N-k)\gamma_k\right].}
$$

$$\rho_k=\gamma_k/\gamma_0$$라고 하면

$$
\operatorname{Var}(\mu)
=\frac{\gamma_0}{N}\left[1+2\sum_{k=1}^{N-1}
\left(1-\frac kN\right)\rho_k\right].
$$

양의 자기상관이 지배적이면 단순 i.i.d. 식보다 평균의 분산이 커질 수 있다. 적절한 수렴 조건에서 장기분산은

$$
\Omega=\gamma_0+2\sum_{k=1}^{\infty}\gamma_k
$$

이며 큰 $$N$$에서 $$\operatorname{Var}(\mu)\approx\Omega/N$$이다. 실제 추정에서는 적절한 시차 절단과 가중을 사용한 HAC 추정 등으로 연결할 수 있다. 이 근사를 사용하려면 의존성에 관한 추가 조건이 필요하다.

같은 전개로 합산 수익률의 분산은

$$
\operatorname{Var}\left(\sum_{i=1}^{N}r_i\right)
=N\gamma_0+2\sum_{k=1}^{N-1}(N-k)\gamma_k.
$$

따라서 자기상관이 있으면 변동성을 무조건 $$\sqrt N$$배 하는 기간 환산도 정확하지 않다.

### 부록 G. 정규표본의 분산 추정은 얼마나 정확한가? — 추가 보충
{: #appendix-g }

표본이 정규 i.i.d.인 경우, 다음 카이제곱 분포 정리를 사용할 수 있다.

$$
Q=\frac{(N-1)s^2}{\sigma_{\mathrm{pop}}^2}\sim\chi^2_{N-1}.
$$

정리 자체의 선형대수적 증명은 여기서 생략하고, 이를 이용한 분산 계산을 전개한다. 자유도 $$\nu$$인 카이제곱변수는 $$\operatorname{Var}(Q)=2\nu$$이므로

$$
\begin{aligned}
\operatorname{Var}(s^2)
&=\left(\frac{\sigma_{\mathrm{pop}}^2}{N-1}\right)^2\operatorname{Var}(Q)\\
&=\frac{\sigma_{\mathrm{pop}}^4}{(N-1)^2}2(N-1)\\
&=\frac{2\sigma_{\mathrm{pop}}^4}{N-1}.
\end{aligned}
$$

따라서

$$
\operatorname{SE}(s^2)=\sigma_{\mathrm{pop}}^2\sqrt{\frac2{N-1}},
\qquad
\widehat{\operatorname{SE}}(s^2)=s^2\sqrt{\frac2{N-1}}.
$$

이 정확한 분산식은 **정규 i.i.d. 가정**에 의존한다. 두꺼운 꼬리나 시간적 의존성이 있으면 그대로 적용할 수 없다. 평균뿐 아니라 위험 추정치 자체에도 불확실성이 있다는 점이 핵심이다.

