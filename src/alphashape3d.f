C     FORTRAN SUBROUTINES IN PACKAGE alphashape3d
C
C     Thomas Lafarge <thowas.lafarge@gmail.com>
C     Beatriz Pateiro-López <beatriz.pateiro@usc.es>
C ==========================================================


C     SUBROUTINE SORT(N,A)
C     SORTS THE INTEGER ARRAY A(I), I=1,2,...,N IN ASCENDING ORDER
C     ALGORTIHM INSERTION SORT (faster for small N)
      SUBROUTINE SORT(N,A)
      IMPLICIT NONE
      INTEGER N,I,J
      INTEGER A(N),X
      DO 30 I=2,N
      X=A(I)
      J=I
   10 J=J-1
      IF(J.EQ.0 .OR. A(J).LE.X) GO TO 20
      A(J+1)=A(J)
      GO TO 10
   20 A(J+1)=X
   30 CONTINUE
      END

C     SUBROUTINE SORTM (ORDM,NR,NC,M)
C     SORTS THE ROWS OF THE INTEGER MATRIX M IN INCREASING ORDER
      SUBROUTINE SORTM (ORDM,NR,NC,M)
      INTEGER NR, NC, I, J
      INTEGER M(NR,NC), ORDM(NR,NC)
      INTEGER A(NC)
      DO I = 1,NR
      A=M(I,:)
      CALL SORT(NC,A)
      DO J=1,NC
      ORDM(I,J)=A(J)
      END DO
      END DO
      RETURN
      END

C     SUBROUTINE Fm123(x,nrx,t1,t2,t3,nt,m123)
C     COMPUTES DETERMINANTS Mijk
      SUBROUTINE Fm123(x,nrx,t1,t2,t3,nt,m123)
      INTEGER nrx, nt
      DOUBLE PRECISION x(nrx,3)
      INTEGER t1(nt), t2(nt), t3(nt)
      DOUBLE PRECISION aux(3,3), m123(nt)
      DOUBLE PRECISION  a1,a2,a3
      INTEGER i
      DO i=1,nt
      aux(1,:)=x(t1(i),:)
      aux(2,:)=x(t2(i),:)
      aux(3,:)=x(t3(i),:)
      a1=aux(1,1)*(aux(2,2)*aux(3,3)-aux(2,3)*aux(3,2))
      a2=-aux(1,2)*(aux(2,1)*aux(3,3)-aux(2,3)*aux(3,1))
      a3=aux(1,3)*(aux(2,1)*aux(3,2)-aux(2,2)*aux(3,1))
      m123(i)=a1+a2+a3
      END DO
      RETURN
      END

C     SUBROUTINE Fmk0(x,nrx,e1,e2,ned,mk0,num)
C     COMPUTES DETERMINANTS Mk0, k=1,2,3 AND NUM.RHO2
      SUBROUTINE Fmk0(x,nrx,e1,e2,ned,mk0,num)
      INTEGER nrx, ned,i
      DOUBLE PRECISION x(nrx,3)
      INTEGER e1(ned),e2(ned)
      DOUBLE PRECISION mk0(ned,3),num(ned)
      DO i=1,ned
      mk0(i,1)=x(e1(i),1)-x(e2(i),1)
      mk0(i,2)=x(e1(i),2)-x(e2(i),2)
      mk0(i,3)=x(e1(i),3)-x(e2(i),3)
      num(i)=mk0(i,1)**2+mk0(i,2)**2+mk0(i,3)**2
      END DO
      RETURN
      END

C     SUBROUTINE Fmij0(x,n,t1,t2,t3,nt,ta,ie,ned,m0,m23,m13,m12,nu2,nu3)
C     COMPUTES DETERMINANTS Mij0, AND NUM.RHO3
      SUBROUTINE Fmij0(x,n,t1,t2,t3,nt,ta,ie,ned,m0,m23,m13,m12,nu2,nu3)
      INTEGER n, nt, ned,i
      DOUBLE PRECISION x(n,3)
      INTEGER t1(nt), t2(nt), t3(nt), ta(nt,3),ie(3*nt)
      DOUBLE PRECISION m0(ned,3),m23(nt),m13(nt),m12(nt)
      DOUBLE PRECISION nu2(ned),nu3(nt)
      DOUBLE PRECISION au1,au2,au3
      DO i=1,nt
      au1=x(t1(i),2)*m0(ie(ta(i,3)),3)
      au2=-x(t2(i),2)*m0(ie(ta(i,2)),3)
      au3=x(t3(i),2)*m0(ie(ta(i,1)),3)
      m23(i)=au1+au2+au3
      au1=x(t1(i),1)*m0(ie(ta(i,3)),3)
      au2=-x(t2(i),1)*m0(ie(ta(i,2)),3)
      au3=x(t3(i),1)*m0(ie(ta(i,1)),3)
      m13(i)=au1+au2+au3
      au1=x(t1(i),1)*m0(ie(ta(i,3)),2)
      au2=-x(t2(i),1)*m0(ie(ta(i,2)),2)
      au3=x(t3(i),1)*m0(ie(ta(i,1)),2)
      m12(i)=au1+au2+au3
      nu3(i)=nu2(ie(ta(i,1)))*nu2(ie(ta(i,2)))*nu2(ie(ta(i,3)))
      END DO
      RETURN
      END

C     SUBROUTINE Fmijk0(x,n,tc,ntc,tca,it,nt,m23,m13,m12,m234,m134,m124)
C     COMPUTES DETERMINANTS Mijk0
      SUBROUTINE Fmijk0(x,n,tc,ntc,tca,it,nt,m23,m13,m12,m234,m134,m124)
      INTEGER n,ntc,nt,i
      DOUBLE PRECISION x(n)
      INTEGER tc(ntc,4),tca(ntc,4),it(4*ntc)
      DOUBLE PRECISION m23(nt),m13(nt),m12(nt)
      DOUBLE PRECISION m234(ntc),m134(ntc),m124(ntc),m123(ntc)
      DOUBLE PRECISION au1,au2,au3,au4
      DO i=1,ntc
      au1=x(tc(i,1))*m23(it(tca(i,4)))
      au2=-x(tc(i,2))*m23(it(tca(i,3)))
      au3=x(tc(i,3))*m23(it(tca(i,2)))
      au4=-x(tc(i,4))*m23(it(tca(i,1)))
      m234(i)=au1+au2+au3+au4
      au1=x(tc(i,1))*m13(it(tca(i,4)))
      au2=-x(tc(i,2))*m13(it(tca(i,3)))
      au3=x(tc(i,3))*m13(it(tca(i,2)))
      au4=-x(tc(i,4))*m13(it(tca(i,1)))
      m134(i)=au1+au2+au3+au4
      au1=x(tc(i,1))*m12(it(tca(i,4)))
      au2=-x(tc(i,2))*m12(it(tca(i,3)))
      au3=x(tc(i,3))*m12(it(tca(i,2)))
      au4=-x(tc(i,4))*m12(it(tca(i,1)))
      m124(i)=au1+au2+au3+au4
      END DO
      RETURN
      END



C     SUBROUTINE int3(ntri2,rf1,rf2,l3,u3)
C     COMPUTES INTERVALS OF ALPHA VALUES FOR TRIANGLES
      SUBROUTINE int3(ntri2,rf1,rf2,l3,u3)
      INTEGER ntri2
      DOUBLE PRECISION rf1(ntri2), rf2(ntri2),l3(ntri2), u3(ntri2)
      INTEGER i
      DO i=1,ntri2
      IF(rf1(i).GT.rf2(i)) THEN
      u3(i)=rf1(i)
      l3(i)=rf2(i)
      ELSE
      u3(i)=rf2(i)
      l3(i)=rf1(i)
      END IF
      END DO
      END

C     SUBROUTINE int2(dup,ned,ntri,itro,l3,u3,l2,u2,tra,aux,ie,nrh,m,ia)
C     COMPUTES INTERVALS OF ALPHA VALUES FOR EDGES
      SUBROUTINE int2(dup,ned,ntri,itro,l3,u3,l2,u2,tra,aux,ie,nrh,m,ia)
      INTEGER ned,ntri
      DOUBLE PRECISION l3(ntri),u3(ntri)
      DOUBLE PRECISION l2(ned),u2(ned),nrh(ned)
      DOUBLE PRECISION m(ned,3)
      DOUBLE PRECISION aa1,aa2,aa3
      INTEGER dup(ned)
      INTEGER itro(3*ntri),tra(ntri,3),aux(3*ntri),ie(3*ntri)
      INTEGER i,j,k,h
      INTEGER ed(3),edo(3),ia(ned)
      k=0
      DO i=1,ned
      ia(i)=0
      l2(i)=1.0E16
      u2(i)=-1.0E16
      DO j=1,dup(i)
      k=k+1
      l2(i)=min(l2(i),l3(itro(k)))
      u2(i)=max(u2(i),u3(itro(k)))
      IF (ia(i).EQ.0) THEN
      DO h=1,3
      ed(h)=tra(itro(k),h)
      edo(h)=ie(ed(h))
      END DO
      IF (aux(k).EQ.1) THEN
      aa1=(m(edo(2),1)+m(edo(3),1))**2
      aa2=(m(edo(2),2)+m(edo(3),2))**2
      aa3=(m(edo(2),3)+m(edo(3),3))**2
      at=nrh(i)-aa1-aa2-aa3
      ELSE IF (aux(k).EQ.2) THEN
      aa1=(m(edo(1),1)-m(edo(3),1))**2
      aa2=(m(edo(1),2)-m(edo(3),2))**2
      aa3=(m(edo(1),3)-m(edo(3),3))**2
      at=nrh(i)-aa1-aa2-aa3
      ELSE IF (aux(k).EQ.3) THEN
      aa1=(-m(edo(1),1)-m(edo(2),1))**2
      aa2=(-m(edo(1),2)-m(edo(2),2))**2
      aa3=(-m(edo(1),3)-m(edo(2),3))**2
      at=nrh(i)-aa1-aa2-aa3
      END IF
      IF (at.GT.0) THEN
      ia(i)=1
      END IF
      END IF
      END DO
      END DO
      END

C     SUBROUTINE int1(dup,nvt,ned,iedo,l2,u2,l1,u1)
C     COMPUTES INTERVALS OF ALPHA VALUES FOR VERTICES
      SUBROUTINE int1(dup,nvt,ned,iedo,l2,u2,l1,u1)
      INTEGER ned,ntri
      DOUBLE PRECISION  l2(ned),u2(ned)
      DOUBLE PRECISION  l1(nvt),u1(nvt)
      INTEGER dup(nvt)
      INTEGER iedo(2*ned)
      INTEGER i,j,k
      k=0
      DO i=1,nvt
      l1(i)=1.0E16
      u1(i)=-1.0E16
      DO j=1,dup(i)
      k=k+1
      l1(i)=min(l1(i),l2(iedo(k)))
      u1(i)=max(u1(i),u2(iedo(k)))
      END DO
      END DO
      END

C     SUBROUTINE edgeSelect(n, x, ed, nb, alpha, nbfinal)
C     SELECTION OF EDGES
      SUBROUTINE edgeSelect(n, x, ed, nb, alpha, nbfinal)
      INTEGER i
      INTEGER n
      INTEGER nbfinal
      DOUBLE PRECISION x(n*3)
      INTEGER nb
      INTEGER ed(nb*2)
      DOUBLE PRECISION dist
      DOUBLE PRECISION alpha
      nbfinal=0
      DO i = 1, nb
      dist = (x(ed(i))- x(ed(i+nb)))**2
      dist = dist+(x(ed(i)+n)-(x(ed(i+nb)+n)))**2
      dist = sqrt(dist+(x(ed(i)+2*n)-(x(ed(i+nb)+2*n)))**2)
      IF (dist < 2*alpha) THEN
      nbfinal=nbfinal+1
      ed(nbfinal)=i
      END IF
      END DO
      END
