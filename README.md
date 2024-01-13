# nTomo.jl

**nTomo.jl** 은 중성자 CT(computer tomography) 의 image reconstruction 을 위한 패키지이다. 

</br>

## 설치

Julia REPL 에서 `]` 를 통해 `Pkg` API 로 진입한 후

```txt
(@v1.10) pkg> add https://github.com/Julia-KAERI/nTomo.jl.git
```

를 통해 설치 할 수 있다.

</br>

## 현재의 기능

- [하나로](https://www.kaeri.re.kr/hanaroe) 의 중성자 영상 장치 (NRF, ENF) 에서 측정한 토모그래피 데이터를 읽는다.
- 데이터의 중간값 필터링 및 선택적 중간값 필터링.
- Filtered back projection(FBP) 를 이용한 3차원 이미지 재구성
- 재구성된 3d volume 을 "tif" image stack 이나 [nrrd](https://en.wikipedia.org/wiki/Nrrd) 포맷으로 저장 할 수 있다.

</br>
