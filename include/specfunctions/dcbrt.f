! cubic root function in double precision
!
      function dcbrt(x)
      implicit real*8 (a - h, o - z)
      dimension c(0 : 23)
      parameter (
     &    p2pow16 = 65536.0d0, 
     &    p2pow48 = 281474976710656.0d0)
      parameter (
     &    p2powm16 = 1 / p2pow16, 
     &    p2powm48 = 1 / p2pow48)
      data c / 
     &    1.5319394088521d-3, -1.8843445653409d-2, 
     &    1.0170534986000d-1, -3.1702448761286d-1, 
     &    6.3520892642253d-1, -8.8106985991189d-1, 
     &    1.0517503764540d0, 4.2674123235580d-1, 
     &    1.5079083659190d-5, -3.7095709111375d-4, 
     &    4.0043972242353d-3, -2.4964114079723d-2, 
     &    1.0003913718511d-1, -2.7751961573273d-1, 
     &    6.6256121926465d-1, 5.3766026150315d-1, 
     &    1.4842542902609d-7, -7.3027601203435d-6, 
     &    1.5766326109233d-4, -1.9658008013138d-3, 
     &    1.5755176844105d-2, -8.7413201405100d-2, 
     &    4.1738741349777d-1, 6.7740948115980d-1 / 
      if (x .eq. 0) then
          dcbrt = 0
          return
      end if
      if (x .gt. 0) then
          w = x
          y = 0.5d0
      else
          w = -x
          y = -0.5d0
      end if
      if (w .gt. 8) then
          do while (w .gt. p2pow48)
              w = w * p2powm48
              y = y * p2pow16
          end do
          do while (w .gt. 8)
              w = w * 0.125d0
              y = y * 2
          end do
      else if (w .lt. 1) then
          do while (w .lt. p2powm48)
              w = w * p2pow48
              y = y * p2powm16
          end do
          do while (w .lt. 1)
              w = w * 8
              y = y * 0.5d0
          end do
      end if
      if (w .lt. 2) then
          k = 0
      else if (w .lt. 4) then
          k = 8
      else
          k = 16
      end if
      u = ((((((c(k) * w + c(k + 1)) * w + 
     &    c(k + 2)) * w + c(k + 3)) * w + 
     &    c(k + 4)) * w + c(k + 5)) * w + 
     &    c(k + 6)) * w + c(k + 7)
      dcbrt = y * (u + 3 * u * w / (w + 2 * u * u * u))
      end
!
