!+++++++++++++++++++++++++++++
! 水平2次元移流分散解析（差分法）
!+++++++++++++++++++++++++++++
program reidai126
      implicit none
      double precision :: c(1510,610,4),cold(1510,610,4)
      double precision :: y(4),fo(4),cs(4)
      double precision :: xmax,ymax,dx,dy,dx2,dy2,vel,al,at,tau,dm,dd1,dd2
      double precision :: shaba,xob1,xob2,tend,dt,time1,o,p,q,v,w,r1,r2,s,u
      integer :: imin,imax,jmin,jmax,jkaz,jc,jj1,jj2,ig1,ig2,nmax,jout1
      integer :: kmin,kmax,jout,k,i,j,n

! 出力ファイルの設定
      open(11,file='output1.csv')
      open(12,file='output2.csv')
      open(13,file='output3.csv')
      open(14,file='output4.csv')
      open(15,file='output5.csv')
      open(16,file='output6.csv')

! 解析に用いる諸数値の設定
      xmax = 1000.0d0                           ! 解析領域(x方向)(m)
      ymax = 500.0d0                            ! 解析領域(y方向)(m)
      dx = 1.0d0                                ! 差分格子間隔(x方向)(m)
      dy = 1.0d0                                ! 差分格子間隔(y方向)(m)
      dx2 = dx**2                               ! 差分格子間隔(x方向)の二乗
      dy2 = dy**2                               ! 差分格子間隔(y方向)の二乗
      imin = 1                                  ! 上流境界の格子点番号
      imax = xmax / dx + 1                      ! 下流境界の格子点番号
      jmin = 1                                  ! 右側境界の格子点番号
      jmax = ymax / dy + 1                      ! 左側境界の格子点番号

      vel = 1.0d-1                              ! 実流速(m/d)
      al = 10.0d0                               ! 縦分散長(m)
      at = al / 10.0d0                          ! 横分散長(m)
      tau =1.0d0                                ! 屈曲率(-)
      dm = 8.64d-5                              ! 分子拡散係数(m^2/d)
      dd1 = al * abs(vel) + dm * tau            ! 分散係数(x方向)(m^2/d)
      dd2 = at * abs(vel) + dm * tau            ! 分散係数(y方向)(m^2/d)

      cs(1) = 1.0d0                             ! PCEの汚染源濃度(mg/L)
      cs(2) = 0.0d0                             ! TCEの汚染源濃度(mg/L)
      cs(3) = 0.0d0                             ! DCEの汚染源濃度(mg/L)
      cs(4) = 0.0d0                             ! CEの汚染源濃度(mg/L)

      shaba = 10.0d0                            ! 汚染源幅(m)
      jkaz = int(shaba / dy)                    ! 汚染源幅分の格子点数
      jc = int(jmax / 2) + 1                    ! 汚染源中心の格子点番号
      jj1 = jc - jkaz / 2                       ! 汚染源中心から右側の汚染源端の格子点番号
      jj2 = jc + jkaz / 2                       ! 汚染源中心から左側の汚染源幅の格子点番号

      xob1 = 100.0d0                            ! 濃度観測地点1のx座標(m)
      xob2 = 400.0d0                            ! 濃度観測地点2のx座標(m)
      ig1 = xob1 / dx + 1                       ! 濃度観測地点1の格子点番号
      ig2 = xob2 / dx + 1                       ! 濃度観測地点2の格子点番号

      tend = 365 * 20                           ! 計算時間(d)
      dt = 5.0d-2                               ! 差分時間間隔(d)
      nmax = tend / dt                          ! 総時間ステップ数

      jout1 = 5 / dt                            ! 濃度観測地点の濃度をファイルに出力する時間ステップ数

      kmin = 1                                  ! 対象物質数の最小値
      kmax = 4                                  ! 対象物質数の最大値（対象物質数）
      fo(1) = 1.0d-3                            ! PCEの一次反応速度定数(1/d)
      fo(2) = 1.0d-3                            ! TCEの一次反応速度定数(1/d)
      fo(3) = 1.0d-4                            ! DCEの一次反応速度定数(1/d)
      fo(4) = 1.0d-4                            ! CEの一次反応速度定数(1/d)
      y(1) = 0.0d0                              ! PCEの算出率(PCEは親物質で生成しないのでゼロ)
      y(2) = 131.39d0 / 165.83d0                ! TCEの算出率(=TCEの分子量/PCEの分子量)
      y(3) = 96.95d0 / 131.39d0                 ! DCEの算出率(=DCEの分子量/TCEの分子量)
      y(4) = 62.5d0 / 96.95d0                   ! CEの算出率(=CEの分子量/DCEの分子量)

! 初期条件
      time1 = 0.0d0
      jout = 0
      do k = imin, kmax
      do i = imin, imax
      do j = jmin, jmax
      c(i,j,k) = 0.0
      if((i == imin).and.(j >= jj1).and.(j <=jj2)) c(i,j,k) = cs(k)
      cold(i,j,k) = c(i,j,k)
      end do
      end do
      end do

! 濃度観測地点の初期時間と初期濃度のファイル出力
      write(11,*) 'elapsed time(days),PCE(mg/L),TCE(mg/L),DCE(mg/L),CE(mg/L)'
      write(12,*) 'elapsed time(days),PCE(mg/L),TCE(mg/L),DCE(mg/L),CE(mg/L)'
      write(11,118) time1,c(ig1,jc,1),c(ig1,jc,2),c(ig1,jc,3),c(ig1,jc,4)
      write(12,118) time1,c(ig2,jc,1),c(ig2,jc,2),c(ig2,jc,3),c(ig2,jc,4)
 118  format(f7.1,',',f13.9,',',f13.9,',',f13.9,',',f13.9)

!///// 計算開始 /////
      do n = 1, nmax                            ! 時間ステップの繰り返し
      time1 = time1 + dt
      jout = jout + 1
       if(jout == jout1) then
       write(*,331) time1
 331   format(f10.2,2x,'days')
       else
       end if

      do k = kmin, kmax                         ! 物質の繰り返し
       do i = imin, imax                        ! x方向の格子点の繰り返し
        do j = jmin, jmax                       ! y方向の格子点の繰り返し

! 上流境界条件 
        if((i == imin).and.(j >= jj1).and.(j <= jj2)) then
        cold(i,j,k) = cs(k)
        go to 113
        else if((i == imin).and.(j < jj1)) then
        cold(i,j,k) = 0.0d0
        go to 113
        else if ((i == imin).and.(j > jj2)) then
        cold(i,j,k) = 0.0d0
        go to 113
        else
        end if
! 下流境界条件
        if(i == imax) then
        cold(i,j,k) = 0.0d0
        go to 113
        else
        end if
! 右側境界条件
        if(j == jmin) then
        cold(i,j-1,k) = cold(i,j+1,k)
        else
        end if
! 左側境界条件
        if(j == jmax) then
        cold(i,j+1,k) = cold(i,j-1,k)
        else
        end if

! 境界以外の濃度の計算
        o = cold(i+1,j,k)
        p = cold(i, j, k)
        q = cold(i-1,j,k)
        v = cold(i,j+1,k)
        w = cold(i,j-1,k)
        r1 = dd1 * (o - 2 * p + q) / dx2
        r2 = dd2 * (v - 2 * p + w) / dy2
        s = vel * (o - q) / 2 / dx
        u = y(k) * fo(k-1) * cold(i,j,k-1) - fo(k) * cold(i,j,k)
        c(i,j,k) = p + (r1 + r2 - s + u) * dt
 113    end do                                ! 次のy方向の格子点へ
       end do                                 ! 次のx方向の格子点へ

! 既知濃度の更新
        do i = imin, imax
        do j = jmin, jmax
        cold(i,j,k) = c(i,j,k)
        end do
        end do

      end do                                  ! 次の物質へ

! 一定時間ごとの濃度観測地点における時間と濃度のファイル出力
        if(jout == jout1) then
        write(11,118) time1,c(ig1,jc,1),c(ig1,jc,2),c(ig1,jc,3),c(ig1,jc,4)
        write(12,118) time1,c(ig2,jc,1),c(ig2,jc,2),c(ig2,jc,3),c(ig2,jc,4)
        jout=0
        else
        end if

      end do                                  ! 次の時間ステップへ
!///// 計算終了 /////

! 等濃度線図を描くための各格子点濃度のファイル出力
      write(13,*) 'x(m),y(m),PCE(mg/L)'
      write(14,*) 'x(m),y(m),TCE(mg/L)'
      write(15,*) 'x(m),y(m),DCE(mg/L)'
      write(16,*) 'x(m),y(m),CE(mg/L)'
      do i = imin, imax,20
      do j = jmin, jmax,10
      write(13,1111) dx*(i-1),dy*(j-1),c(i,j,1)
      write(14,1111) dx*(i-1),dy*(j-1),c(i,j,2)
      write(15,1111) dx*(i-1),dy*(j-1),c(i,j,3)
      write(16,1111) dx*(i-1),dy*(j-1),c(i,j,4)
 1111 format(f7.1,',',f7.1,',',f13.9)
      end do
      end do

! 出力ファイルを閉じる
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
      close(16)

      stop
end program reidai126