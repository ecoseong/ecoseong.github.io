---
layout: post
title: "[경제수학] Chapter 1. 경제수학이란"
date: 2026-08-05
categories: [경제수학]
tags: [경제수학, 김완진, 정학, 비교정학, 동학, 이윤극대화]
math: true
---

## 1. 경제학 연구에서 최적화와 균형

경제학 연구에서 핵심적인 개념은 **최적화**와 **균형**이다.

어떤 특정한 경제문제를 설명하기 위해 경제학자가 가장 먼저 하는 일은, 그 문제를 설명하는 데 적합한 **경제모형**을 설정하는 것이다.

- **균형분석(정학, statics)**: 경제모형의 균형을 찾고, 그 균형상태에서 어떤 결과가 나타나는지 분석하는 방법
- **비교정학(comparative statics)**: 정책의 시행이나 다른 조건이 변화할 때 나타나는 새로운 균형을, 최초의 균형과 비교하여 분석하는 방법
- **동학(dynamics)**: 그 변화의 과정이 시간적으로 어떻게 나타나는가를 분석하는 방법

## 2. 내생변수와 외생변수

가장 간단한 예로 수요공급모형을 통해 설명해보자. 수요곡선과 공급곡선은 아래 그림처럼 존재하고,

$$
Q_D = D(p, M), \qquad Q_S = S(p, r)
$$

로 표현할 수 있다.

![그림 1 수요공급과 균형](/assets/images/econmath/fig1_supply_demand.png)
*&lt;그림 1&gt; 수요공급과 균형*

이 두 곡선이 만나는 지점이 **균형조건식**이며, 이때의 균형가격, 균형수요량, 균형공급량은 방정식에 의해 결정되는 변수로서 **내생변수(endogenous variable)**가 된다.

반면 소비자의 소득 $M$, 생산요소가격 $r$은 이 모형에서 외생적으로 주어진 값으로, **외생변수(exogenous variable)** 혹은 **파라미터(parameter)**라고 부른다.

위 그래프처럼 모형을 설정하고, 그 모형의 균형을 찾아 경제현상을 설명하는 분석방법이 바로 **정학**이다.

## 3. 비교정학 분석: 수요공급모형 예시

수요함수와 공급함수를 다음과 같이 선형으로 설정하자.

$$
Q_D = M - ap \quad (a > 0)
$$

$$
Q_S = cp - dr \quad (c, d > 0)
$$

균형조건 $Q_D = Q_S$로부터 균형가격은 다음과 같이 결정된다.

$$
p^0 = \frac{M + dr}{a + c}
$$

**비교정학분석**은 $M$, $r$과 같은 외생변수의 변화가 균형값 $p^0$에 미치는 영향을 연구하는 것이다. 이를 위해서는 $p^0$를 $M$, $r$의 함수로 표현하고, 이 함수가 $M$ 혹은 $r$에 대해 증가함수인지 감소함수인지를 판단해야 한다. 여기에 **미분법**이 유용하게 사용된다.

$$
\frac{dp^0}{dM} = \frac{1}{a+c} > 0
$$

$p^0$의 도함수가 양(+)이므로, $p^0$는 $M$의 증가함수임을 알 수 있다. 마찬가지로 $r$에 대해 미분하면 $p^0$가 $r$에 대해 증가함수인지 감소함수인지도 알 수 있다.

## 4. 최적화모형: 이윤극대화

최적화모형의 대표적인 예로 **이윤극대화** 문제가 있다.

이윤 $\pi$는 총수입에서 총비용을 뺀 것으로,

$$
\pi = pQ - C(Q)
$$

로 표현할 수 있다.

![그림 2 이윤극대화](/assets/images/econmath/fig2_profit_max.png)
*이윤극대화: 총수입·총비용곡선, 총이윤곡선, 한계수입·한계비용곡선*

이 $\pi$를 극대화하기 위해서는 미분법을 이용하면 된다. **1계조건(FOC)**과 **2계조건(SOC)**을 통해 이윤이 언제 극대화되는지 구할 수 있는데, 이 과정에서 필요한 것이 바로 미분이다.

이 책은 이러한 경제학 공부에 필요한 수학적 내용을 소개하는 것을 목적으로 한다.

---

## 부록: 그림 재현 Python 코드

위 두 그림은 아래 코드로 그린 것이다 (원본 교재 그림의 개형을 참고하여 재구성한 것으로, 정확한 수치가 아닌 개형 위주의 근사임).

### 그림 1. 수요공급과 균형

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# 한글 폰트 등록 (환경에 맞게 경로/이름 교체: macOS는 'AppleGothic', Windows는 'Malgun Gothic')
fm.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

fig, ax = plt.subplots(figsize=(5, 5))

q = np.linspace(0, 10, 100)
Pd = 8 - 0.6 * q   # 수요곡선 (우하향)
Ps = 1 + 0.6 * q   # 공급곡선 (우상향)

idx = np.argmin(np.abs(Pd - Ps))
Q0, P0 = q[idx], Pd[idx]
P1, P2 = P0 + 1.5, P0 - 1.5
Qa, Qc = Q0 - 1.5, Q0 + 1.5

ax.plot(q, Pd, color="black")
ax.plot(q, Ps, color="black")

for p in [P0, P1, P2]:
    ax.hlines(p, 0, 10, linestyles="dotted", color="gray", linewidth=0.8)
for qx in [Qa, Q0, Qc]:
    ax.vlines(qx, 0, 8.5, linestyles="dotted", color="gray", linewidth=0.8)

ax.text(-0.6, P1, "$P_1$", va="center")
ax.text(-0.6, P0, "$P_0$", va="center")
ax.text(-0.6, P2, "$P_2$", va="center")
ax.text(Qa, -0.5, "a", ha="center")
ax.text(Q0, -0.5, "b", ha="center")
ax.text(Qc, -0.5, "c", ha="center")
ax.text(1.2, 7.5, "수요", fontsize=11)
ax.text(8.3, 7.5, "공급", fontsize=11)

ax.set_xlim(0, 10)
ax.set_ylim(0, 9)
ax.set_xlabel("수량", loc="right")
ax.set_ylabel("가격", loc="top", rotation=0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

plt.title("<그림 1> 수요공급과 균형", pad=15)
plt.tight_layout()
plt.savefig("fig1_supply_demand.png", dpi=150)
```

### 그림 2. 이윤극대화 (총수입·총비용 / 총이윤 / 한계수입·한계비용)

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

fm.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

Q = np.linspace(0.1, 10, 400)

TC = 1.2 * Q + 0.15 * (Q - 5) ** 3 / 5 + 6      # 총비용곡선 (S자형)
TR = 9 * Q - 0.45 * Q**2                         # 총수입곡선 (오목)
profit = TR - TC

dQ = Q[1] - Q[0]
MC = np.gradient(TC, dQ)
MR = np.gradient(TR, dQ)

fig, axes = plt.subplots(3, 1, figsize=(6, 10), sharex=True)

axes[0].plot(Q, TC, color="black")
axes[0].plot(Q, TR, color="steelblue")
axes[0].set_ylabel("TR, TC\n(총수입, 총비용)")
axes[0].set_title("(i) 총수입, 총비용곡선", fontsize=10)
axes[0].text(Q[-1] + 0.1, TC[-1], "$TC(Q)$", va="center", fontsize=9)
axes[0].text(Q[-1] + 0.1, TR[-1], "$TR(Q)$", va="center", color="steelblue", fontsize=9)

axes[1].axhline(0, color="gray", linewidth=0.8)
axes[1].plot(Q, profit, color="firebrick")
axes[1].set_ylabel("π (이윤)")
axes[1].set_title("(ii) 총이윤곡선", fontsize=10)
axes[1].text(Q[-1] + 0.1, profit[-1], "$\\pi(Q)$", va="center", color="firebrick", fontsize=9)

axes[2].plot(Q, MC, color="black")
axes[2].plot(Q, MR, color="steelblue")
axes[2].set_ylabel("MR, MC\n(한계수입, 한계비용)")
axes[2].set_xlabel("Q (산출량)")
axes[2].set_title("(iii) 한계수입, 한계비용곡선", fontsize=10)
axes[2].text(Q[-1] + 0.1, MC[-1], "$MC(Q)$", va="center", fontsize=9)
axes[2].text(Q[-1] + 0.1, MR[-1], "$MR(Q)$", va="center", color="steelblue", fontsize=9)

for ax in axes:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
axes[0].set_xlim(0, 10.8)

plt.tight_layout()
plt.savefig("fig2_profit_max.png", dpi=150)
