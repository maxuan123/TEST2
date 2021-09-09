!****************************************************************************
!  Covariance.f90                                                           !
!                                                                           !   
!   @2014 North China Electric Power University. All rights reserved.       !
!	Author : Ma Xubo                                                        !
!                                                                           !
!	Date : 2019.5.11                                                        !
!                                                                           !
!****************************************************************************
!                                                                           !
!  MODULE: ACEP                                                             !
!                                                                           !
!  PURPOSE:  Process cross section for MCNP or RMC code.                    !
!                                                                           !
!****************************************************************************
    
MODULE ACEPModule

USE m_real_initial 
USE M_physicsp
USE M_general,             ONLY : ENDFIO,GetAvailableUnit
USE LibData_module,        ONLY : One_line_read1,Total_MF_number
USE LibData_module,        ONLY : ENDFDataTYPE
USE M_Input,               ONLY: InputDataParameters
USE time,                  ONLY: timer

IMPLICIT NONE

!---------------------------------------------------------------------------
! Define MultiGroup parameters
!---------------------------------------------------------------------------
TYPE(timer)               :: tim

TYPE ACEPData_Bank
     REAL(DBL) :: ZA,AWR                                    !ZA=1000*Z+A，AWR原子中子重量比
     INTEGER :: NIS                                         !Number of isotopes in the material    !同位素的数目
     REAL(KIND=DBL) :: ZAI                                  !designation for an isotope
     REAL(KIND=DBL) :: ABN                                  !Abundance of an isotope in the material  ！材料中同位素的丰度
     INTEGER :: LFW                                         !Flag indicating whether average fission widths are given in the unresolved
     
     
CONTAINS
    PROCEDURE :: ACEProcessMain   => ACEProcessMain_sub
    !PROCEDURE :: ACEPData => ACEPData_sub
    !PROCEDURE :: File4DataOutput => File4DataOutput_sub
!PROCEDURE :: File5DataOutput  => File5DataOutput_sub
    PROCEDURE :: File5DataProcess => File5DataProcess_sub
    PROCEDURE :: File6DataProcess => File6DataProcess_sub
    
END TYPE ACEPData_Bank



    CONTAINS
    
 SUBROUTINE ACEProcessMain_sub(Self,LibIso,FN,FN2,Input,NumMix,NumIso)
!-------------------------------------------------------------------------!
! 该程序主要功能是读取endf原始评价数据库以及由purr模块生成的pendf         !
! 数据库分别将file1-file6中的数据经过处理变成ace格式的数据库，供mcnp使用  !
! file1中读取endf，处理的是平均裂变中子数                                 !
! file2中主要读取pendf库，处理不可分辨数据                                !
! file3中读取pendf库，主要提取能量框架，截面，反应道等信息                !
! file4中读取endf，处理不同能量下的角分布，输出成概率表形式               !
! file5中读取endf，处理能量分布，输出成概率表形式                         !
! file6中读取endf，处理能量角度分布                                       !
!-------------------------------------------------------------------------!
IMPLICIT NONE
CLASS(ACEPData_Bank):: Self
TYPE(ENDFDataTYPE),INTENT(IN) :: LibIso
TYPE(InputDataParameters),INTENT(IN) :: Input
CHARACTER(*) ,INTENT(IN):: FN                               
CHARACTER(*) ,INTENT(IN):: FN2                              ! 评价核数据库原始的数据库名字
INTEGER,INTENT(IN) :: NumMix                                ! MIX 的序号
INTEGER,INTENT(IN) :: NumIso                                ! MIX 中核素的序号

! 局部变量定义
TYPE(ENDFDataTYPE) :: LibIsoFile33 
TYPE(ENDFDataTYPE) :: LibIso1
TYPE(ENDFDataTYPE) :: LibIso2
TYPE(ENDFDataTYPE) :: LibIso3
TYPE(ENDFDataTYPE) :: LibIso5  !放file5
TYPE(ENDFDataTYPE) :: LibIso6  !放file6

TYPE(ENDFDataTYPE) :: LibIso4  !放file4
INTEGER :: NP,NE,TotalPoint,TableNum,DataNum,LCT,mt18,total_number
INTEGER ::MT_total_number,MT_total_number1,MT_total_number5,MT_total_number6
INTEGER :: i,j,k,n,a,t,l,ii,nn,kk,jj,flag              !loop index
INTEGER :: ptr,ptr22,ptr3,ptr33,ptr4,ptr44,ptr444,ptr5,ptr6,ptr455,ptr_dlw,ptr_dlw2,ptr_ldlw       !ptr为pointer
INTEGER :: tot,num
INTEGER,DIMENSION(100)::axis
REAL(KIND=DBL)::mev
REAL(KIND=DBL)::shake
INTEGER,DIMENSION(:),ALLOCATABLE::E
INTEGER,DIMENSION(:),ALLOCATABLE::tyr                       !存放tyr数据
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::xss3               !存放file3中截面，能量框架，平均热数，MTR,TYR等数据
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::xss5
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::xss6
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::coeff
INTEGER,DIMENSION(:),ALLOCATABLE::mtt
INTEGER,DIMENSION(:),ALLOCATABLE::mtr
INTEGER,DIMENSION(:),ALLOCATABLE::mtr2
INTEGER::MT,status,mt_num,m5,m6,num2,num3                                   !mt_num去掉弹性散射的反应道数
!INTEGER ::status 
INTEGER,DIMENSION(0:500) :: Loca                            !Loca每一个反应道中每一部分数据起始位置
INTEGER,DIMENSION(0:100) :: Local                           !Local每一个反应道的起始位置
INTEGER,DIMENSION(100) :: SumLoca                           !LAND中的起始位置
INTEGER,DIMENSION(:),ALLOCATABLE :: m                       !计数用
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: pSumL            !低能区概率之和
REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: angle2         !angle2存放倒栈后的角度
REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE :: cprob2
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: energy_point
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: cross_section
REAL(KIND=DBL):: error,x(70),xm,y(24),b,c,aco(7000),cprob(7000),cumm(7000),angle(7000)             !angle临时存放角度
REAL(KIND=DBL)::prob1,prob2,probm,mmm,dm,ym,yt,xl,yl,check,dco,diff
REAL(KIND=DBL)::delt,cx,cxx,fx,ex,et,dy,f,mf6,count,AWR,xxx
INTEGER ::nxs(16),jxs(32),jx,count2


error=0.001;mev=1000000;shake=100000000
mev=1000000;t=1;mt_num=0;
flag=1


 CALL LibIso1%ENDF_read(FN2,1)
 CALL LibIso2%FILE2153Read(FN)
 CALL LibIso3%ENDF_read(FN2,3) 
 CALL LibIso4%ENDF_read(FN2,4)
 CALL LibIso5%ENDF_read(FN2,5)
 CALL LibIso6%ENDF_read(FN2,6)

 
          DO i=1,4     
              IF(3 == LibIso%MF_number(i) ) THEN
                  MT_total_number = LibIso%Total_MT_number(i)  
              ENDIF
          ENDDO
          DO i=1,5     
              IF(5 == LibIso5%MF_number(i) ) THEN
                  MT_total_number5 = LibIso5%Total_MT_number(i)  
              ENDIF
          ENDDO
          DO i=1,6     
              IF(6 == LibIso6%MF_number(i) ) THEN
                  MT_total_number6 = LibIso6%Total_MT_number(i)  
              ENDIF
          ENDDO
            NP=LibIso%FILE3_MT(1)%NP
            ! head imformation
            nxs(2)= LibIso4%FILE4_MT(1)%ZA
           
            
            nxs(7)=0
            DO i=9,15
                nxs(i)=0
            ENDDO
            nxs(16)=0         !0是正常中子产生光子，-1是不产生光子
            jxs(1)=1
            ALLOCATE(xss3(20000000))
            ALLOCATE(energy_point(NP))
            ALLOCATE(cross_section(NP))
  		  DO i=1,MT_total_number
              MT=LibIso%FILE3_MT(i)%MT          
                IF(LibIso%FILE3_MT(i)%MT.eq.1) THEN
                    PRINT*,NP
                    count=0
                    count2=0
                    !DO j=1,NP-count
                    !    IF((abs(LibIso%FILE3_MT(i)%Energy_point(j+1+count)-LibIso%FILE3_MT(i)%Energy_point(j+count))).le.1) THEN
                    !        PRINT*,J
                    !        xss3(j)=LibIso%FILE3_MT(i)%Energy_point(j+count+1)
                    !        count=count+1
                    !    ELSE              
                    !        xss3(j)=LibIso%FILE3_MT(i)%Energy_point(j+count)
                    !    ENDIF
                    !ENDDO
                    DO j=1,NP
                        !IF(LibIso%FILE3_MT(i)%Energy_point(j).EQ.1090000) THEN
                        !    CONTINUE
                        !ELSE
                            count=count+1
                            xss3(count)=LibIso%FILE3_MT(i)%Energy_point(j)/mev 
                            energy_point(count)=LibIso%FILE3_MT(i)%Energy_point(j)/mev 
                        !ENDIF
                    ENDDO
                    nxs(3)=count

                    DO j=1,NP
                        !IF(LibIso%FILE3_MT(i)%Energy_point(j).EQ.1090000) THEN
                        !    CONTINUE
                        !ELSE
                            count2=count2+1
                            xss3(count2+count)=LibIso%FILE3_MT(i)%Cross_section(j)
                            cross_section(count2)=LibIso%FILE3_MT(i)%Cross_section(j)
                        !ENDIF       
                    ENDDO
                ENDIF
          ENDDO
          DO i=1,MT_total_number
              MT=LibIso%FILE3_MT(i)%MT          
                IF(LibIso%FILE3_MT(i)%MT.eq.102) THEN  
                    COUNT2=0
                    DO j=1,NP
                        
                        !IF(LibIso%FILE3_MT(i)%Energy_point(j).EQ.1090000) THEN
                        !    CONTINUE
                        !ELSE
                            count2=count2+1
                            xss3(count2+2*count)=LibIso%FILE3_MT(i)%Cross_section(j)                 !ESZ BLOCK
                        !ENDIF
                    ENDDO
                ELSEIF(LibIso%FILE3_MT(i)%MT.eq.2) THEN
                    COUNT2=0
                    DO j=1,NP
                        
                        !IF(LibIso%FILE3_MT(i)%Energy_point(j).EQ.1090000) THEN
                        !    CONTINUE
                        !ELSE
                            count2=count2+1
                            xss3(count2+3*count)=LibIso%FILE3_MT(i)%Cross_section(j)
                        !ENDIF
                    ENDDO
                ELSEIF(LibIso%FILE3_MT(i)%MT.eq.301) THEN    !heating factors?
                    COUNT2=0
                    DO j=1,NP
                        !IF(LibIso%FILE3_MT(i)%Energy_point(j).EQ.1090000) THEN
                        !    CONTINUE
                        !ELSE
                            count2=count2+1
                            xss3(count2+4*count)=LibIso%FILE3_MT(i)%Cross_section(j)/ &
                                 LibIso%FILE3_MT(1)%Cross_section(j)/mev
                        !ENDIF
                    ENDDO
                ENDIF 
          ENDDO
          ptr=5*count
          
          !---------------------------
          !NU DATA 平均裂变中子数等信息
          !---------------------------
          count2=0
          MT_total_number1=LibIso1%Total_MT_number(1)
          ALLOCATE(LibIso1%FILE1%MT1(MT_total_number1),stat=status)
          DO i=1,MT_total_number1
              IF(LibIso1%FILE1%MT1(i)%MT.eq.456) THEN
                  jxs(2)=ptr+1
                  IF(LibIso1%FILE1%LNU456.eq.2) THEN
                      count2=LibIso1%FILE1%NP456*LibIso1%FILE1%Intflag456+3
                      xss3(ptr+1)=-count2
                      xss3(ptr+2)=LibIso1%FILE1%LNU456
                      xss3(ptr+3)=0
                      xss3(ptr+4)=LibIso1%FILE1%NP456
                      IF(LibIso1%FILE1%NR456.ne.1.or.LibIso1%FILE1%Intflag456.ne.2) THEN
                          PRINT*,'m!=1  jnt!=2'
                          STOP
                      ELSE  
                        DO j=1,LibIso1%FILE1%NP456
                            xss3(ptr+4+j)=LibIso1%FILE1%EintNuPrompt(2*j-1)/mev
                            xss3(ptr+4+LibIso1%FILE1%NP456+j)=LibIso1%FILE1%EintNuPrompt(2*j)
                        ENDDO
                      ENDIF
                      ptr=ptr+4+2*LibIso1%FILE1%NP456
                  ELSEIF(LibIso1%FILE1%LNU456.eq.2) THEN
                      print*,'lnu=2,can not be processed!'
                      STOP
                  ENDIF
              ENDIF
          ENDDO
          DO i=1,MT_total_number1
              IF(LibIso1%FILE1%MT1(i)%MT.eq.452) THEN    
                  
                  xss3(ptr+1)=LibIso1%FILE1%LNU452
                  xss3(ptr+2)=0
                  xss3(ptr+3)=LibIso1%FILE1%NP452
                  IF(LibIso1%FILE1%NR452.ne.1.or.LibIso1%FILE1%Intflag452.ne.2) THEN
                      PRINT*,'m!=1  jnt!=2'
                      STOP
                  ELSE  
                    DO j=1,LibIso1%FILE1%NP452
                        xss3(ptr+3+j)=LibIso1%FILE1%EintNutotal(2*j-1)/mev
                        xss3(ptr+3+LibIso1%FILE1%NP452+j)=LibIso1%FILE1%EintNutotal(2*j)
                    ENDDO
                  ENDIF
                  ptr=ptr+3+2*LibIso1%FILE1%NP452
              ENDIF
          ENDDO    
              
            ALLOCATE(mtt(MT_total_number))
            ALLOCATE(mtr(MT_total_number))
            DO i=1,MT_total_number
                MT=LibIso%FILE3_MT(i)%MT
                IF(MT.eq.1.or.MT.eq.2.or. &
                 MT.eq.19.or.MT.eq.20.or. &
                 MT.eq.21.or.MT.eq.26.or. &
                 MT.eq.27.or.MT.eq.38.or. &
                 MT.eq.101.or.MT.eq.120.or. &                           !MTR BLOCK
                 MT.eq.151.or. &
                 (MT.gt.207.and.MT.lt.444)) THEN
                  CONTINUE
                ELSEIF(MT.eq.4) THEN
                    !满足一下条件时才处理mt=4
                    IF(LibIso2%MF2153%iinel.eq.0.and.LibIso2%MF2153%iabso.gt.0) THEN
                        LibIso2%MF2153%iinel=mod(int(LibIso2%MF2153%iabso),1000)
                    ENDIF
                    IF(LibIso2%MF2153%iinel.eq.4) THEN
                        mt_num=mt_num+1
                        !xss3(mt_num+5*NP)=MT
                        mtt(mt_num)=i
                        mtr(mt_num)=MT
                    ELSE
                        CONTINUE
                    ENDIF
                ELSE
                    mt_num=mt_num+1
                    !xss3(mt_num+5*NP)=MT
                    mtt(mt_num)=i
                    mtr(mt_num)=MT
                  ENDIF     
            ENDDO
            DO i=1,mt_num
                ii=mtt(i)
                IF(mtr(i).eq.4) THEN
                    DO j=i+1,mt_num
                        mtr(j-1)=mtr(j)      !如果出现mt=4，则把4放在最后一个位置
                        mtt(j-1)=mtt(j)
                    ENDDO
                    mtr(mt_num)=4
                    mtt(mt_num)=ii
                    
                    EXIT
                ELSE
                    CONTINUE
                ENDIF
            ENDDO
            nxs(4)=mt_num
            jxs(3)=ptr+1
            jxs(4)=ptr+mt_num+1
            jxs(5)=ptr+2*mt_num+1
            jxs(6)=ptr+3*mt_num+1
            jxs(7)=ptr+4*mt_num+1
            DO i=1,mt_num
                xss3(ptr+i)=mtr(i)
                xss3(ptr+mt_num+i)=LibIso%FILE3_MT(mtt(i))%QI/mev
            ENDDO
            ALLOCATE(tyr(mt_num),stat=status)
            DO j=1,mt_num                                                                           !TYR Block
                MT=xss3(j+ptr)
                IF(MT==5)THEN
                    tyr(j)=101
                ELSE IF(MT==16)THEN                                                                 
                    tyr(j)=2
                ELSE IF(MT==17)THEN
                    tyr(j)=3
                ELSE IF((MT==18).or.(MT==19).or.(MT==20).or.(MT==21).or.(MT==38))THEN
                    tyr(j)=19
                ELSE IF(MT==37) THEN
                    tyr(j)=4
                ELSE IF((MT>=51).and.(MT<=90))THEN
                    tyr(j)=1
                ELSE IF((MT==91).or.(MT==22).or.(MT==28).or.(MT==32).or.(MT==33))THEN
                    tyr(j)=1
                ELSE IF((MT>100).or.(MT.eq.4)) THEN
                    tyr(j)=0
                ENDIF
                IF(MT.ne.18) THEN
                    xss3(ptr+2*mt_num+j)=-tyr(j)     !初步取负，实际上取决去LCT,质心系取负
                ELSE
                    xss3(ptr+2*mt_num+j)=tyr(j)
                ENDIF
            ENDDO
            DEALLOCATE(tyr)
            ptr3=ptr+3*mt_num
            ptr33=ptr3
            ptr3=ptr3+mt_num     !位置指针
            
            xss3(ptr33+1)=1
            DO k=1,mt_num
                IF(LibIso%FILE3_MT(mtt(k))%MT.eq.18) THEN
                    jxs(21)=ptr3+1
                ENDIF
                DO j=1,LibIso%FILE3_MT(1)%NP
                    IF(LibIso%FILE3_MT(mtt(k))%Energy_point(1)/mev.eq.energy_point(j)) THEN
                        !IF(j.ne.1) THEN
                        !    xss3(ptr3+1)=j      !每一个反应道第一个能量点对应完整能量框架的第几个能量点 
                        !ELSE
                            xss3(ptr3+1)=j                         
                        !ENDIF
                        EXIT
                    ELSE
                        CONTINUE
                    ENDIF
                ENDDO
                IF(xss3(ptr3+1).EQ.1) THEN
                    xss3(ptr3+2)=count
                    count2=0
                    DO j=1,LibIso%FILE3_MT(mtt(k))%NP
                        !IF(LibIso%FILE3_MT(mtt(k))%Energy_point(j).eq.1090000) THEN
                        !    CONTINUE
                        !ELSE
                            count2=count2+1
                            IF(LibIso%FILE3_MT(mtt(k))%MT.eq.444) THEN
                                xss3(ptr3+2+count2)=LibIso%FILE3_MT(mtt(k))%Cross_section(j)/mev
                            ELSE
                                xss3(ptr3+2+count2)=LibIso%FILE3_MT(mtt(k))%Cross_section(j)
                            ENDIF
                        !ENDIF
                    ENDDO
                    ptr3=ptr3+2+count2
                    IF(k.ge.2) THEN
                        !IF(LibIso%FILE3_MT(mtt(k-1))%NP.EQ.count+1) THEN
                        !    xss3(ptr33+k)=xss3(ptr33+k-1)+count+2
                        !ELSE
                        !    xss3(ptr33+k)=xss3(ptr33+k-1)+LibIso%FILE3_MT(mtt(k-1))%NP+2    !每一个反应道的位置指针
                        !ENDIF
                        IF(k.lt.mt_num) THEN
                            xss3(ptr33+k+1)=xss3(ptr33+k)+count2+2
                        ENDIF
                    ENDIF
                ELSE
                    count2=0
                    xss3(ptr3+2)=LibIso%FILE3_MT(mtt(k))%NP
                    DO j=1,LibIso%FILE3_MT(mtt(k))%NP
                        !IF(LibIso%FILE3_MT(mtt(k))%Energy_point(j).eq.1090000) THEN
                        !    CONTINUE
                        !ELSE
                            count2=count2+1
                            IF(LibIso%FILE3_MT(mtt(k))%MT.eq.444) THEN     !mt=444时，数据除以1000000
                                xss3(ptr3+2+count2)=LibIso%FILE3_MT(mtt(k))%Cross_section(j)/mev
                            ELSE
                                xss3(ptr3+2+count2)=LibIso%FILE3_MT(mtt(k))%Cross_section(j)
                            ENDIF
                        !ENDIF
                    ENDDO
                    xss3(ptr3+2)=count2
                    ptr3=ptr3+2+count2
                    IF(k.lt.mt_num) THEN
                        xss3(ptr33+k+1)=xss3(ptr33+k)+count2+2    !每一个反应道的位置指针
                    ENDIF
                ENDIF
            ENDDO
            DEALLOCATE(mtt)
 ! WRITE(ENDFIO%OUTP,'(4ES21.11)')(xss3(i),i=1,ptr3)
  
  ! 处理File4 的数据
  
          
          DO i=1,5      
              IF(4 == LibIso4%MF_number(i) ) THEN
                  MT_total_number = LibIso4%Total_MT_number(i) 
              ENDIF
          ENDDO
          ptr4=ptr3     !指针用改用ptr4表示
            total_number=1
            ptr44=ptr4
            ALLOCATE(mtr2(100))
            mtr2(1)=2
            j=1
            DO i=1,mt_num
                IF(mtr(i).ge.5.and.mtr(i).le.91) THEN
                    total_number=total_number+1
                    j=j+1
                    mtr2(j)=mtr(i)
                ENDIF
            ENDDO
            jxs(8)=ptr44+1
            xss3(ptr44+1)=1
            DO i=2,total_number
                DO j=1,MT_total_number6
                    IF(mtr2(i).eq.LibIso6%FILE6_MT(j)%MT) THEN
                        xss3(ptr44+i)=-1
                    ENDIF
                    
                ENDDO
                IF(mtr2(i).eq.18) THEN
                    xss3(ptr44+i)=0
                ENDIF
            ENDDO
            IF(LibIso4%FILE4_MT(1)%MT.eq.2) THEN
                nxs(5)=total_number-1
            ELSE
                nxs(5)=total_number
            ENDIF
            ptr44=ptr4
           ptr4=ptr4+total_number
           
           ptr444=ptr4
           !------------------------------------!
           !            AND Block               !
           !------------------------------------!
           count=0
           
            jx=0
           DO n=1,MT_total_number

               DO jj=1,total_number
                   IF((LibIso4%FILE4_MT(n)%MT.ne.mtr2(jj)).and.(LibIso4%FILE4_MT(n)%MT.ne.2)) THEN
                       CONTINUE
                   ELSE
                       mf6=0
                       DO nn=1,MT_total_number6
                           IF(LibIso4%FILE4_MT(n)%MT.eq.LibIso6%FILE6_MT(nn)%MT) THEN
                               mf6=1                !mf=6出现的反应道留给mf6处理
                               exit                                                                                       
                           ELSE
                               mf6=0
                           ENDIF
                       ENDDO
                       IF(mf6.ne.0) THEN
                           CONTINUE
                       ELSE
                            LCT=LibIso4%FILE4_MT(n)%LCT
                            NE = LibIso4%FILE4_MT(n)%NE            !NE个需要勒让德转化成概率表
                            TableNum=LibIso4%FILE4_MT(n)%TableNum  !TableNum个概率比，无需再转化
                            IF(NE+TableNum.eq.0) THEN
                                CONTINUE
                            ELSE
                                !ALLOCATE(Loca(NE+TableNum+1),stat=status)
                                ALLOCATE(m(NE+1),stat=status)
                                
                                ALLOCATE(angle2(NE+1000,400),stat=status)
                                ALLOCATE(pSumL(NE+1000),stat=status)
                                if(LibIso4%FILE4_MT(n)%MT.ne.2) then
                                    Loca(1)=2*(NE+TableNum)+1+count
                                else
                                    Loca(1)=2*(NE+TableNum)+2+count
                                endif
                                print*,'ne=',ne
                                DO j=1,NE
                                     TotalPoint = LibIso4%FILE4_MT(n)%E(j)%TotalPoint              !TotalPoint阶数也是系数的个数
                                     !倒栈法求角度
                                     x(1)=1;x(2)=-1;c=x(1);b=x(2);k=2;m(j)=1;angle(1)=-1
                                  04  continue
                                      xm=(b+c)/2
                                      if(abs(func(xm,LibIso4%FILE4_MT(n)%E(j)%Coeffi_a,TotalPoint)- &
                                         (func(c,LibIso4%FILE4_MT(n)%E(j)%Coeffi_a,TotalPoint)+ &
                                         func(b,LibIso4%FILE4_MT(n)%E(j)%Coeffi_a,TotalPoint))/2)<error) go to 05
                                      k=k+1
                                      x(k)=x(k-1)
                                      x(k-1)=xm
                                      b=x(k)
                                      c=x(k-1)
                                      go to 04
                                  05  continue
                                      if(k==2) go to 06
                                      m(j)=m(j)+1
                                      angle(m(j))=c
                                      k=k-1
                                      b=c
                                      c=x(k-1)
                                      go to 04
                                  06  continue
                                      m(j)=m(j)+1
                                      angle(m(j))=1
                                      IF(m(j).eq.2) THEN
                                          Loca(j+1)=Loca(j)+3*3+2
                                      ELSE
                                          Loca(j+1)=Loca(j)+3*m(j)+2
                                      ENDIF
                                      do i=1,m(j)
                                          angle2(i,j)=angle(i)
                                      enddo
                                ENDDO
                                DO j=1,TableNum
                                    DataNum=LibIso4%FILE4_MT(n)%Table(j)%DataNum
                                    Loca(j+1+NE)=Loca(j+NE)+3*DataNum+2
                                ENDDO
                                IF(jx.eq.0) THEN
                                    jxs(9)=ptr4+1
                                    jx=1
                                ENDIF
                                xss3(ptr4+1)=NE+TableNum
                                DO j=1,NE
                                    xss3(ptr4+1+j)=LibIso4%FILE4_MT(n)%E(j)%E/mev
                                ENDDO
                                DO j=1,TableNum
                                    xss3(ptr4+1+NE+j)=LibIso4%FILE4_MT(n)%Table(j)%E/mev
                                ENDDO
                                 !输出每一个概率表的起始位置
                                IF(LCT.eq.1) THEN 
                                    DO j=1,NE+TableNum
                                        xss3(ptr4+1+NE+TableNum+j)=Loca(j)    !32 equiprobable bin distribution
                                    ENDDO
                                ELSEIF(LCT.eq.2) THEN
                                    DO j=1,NE+TableNum
                                        xss3(ptr4+1+NE+TableNum+j)=-Loca(j)   !tabulated angular distribution
                                    ENDDO
                                ENDIF
                                ptr4=ptr4+1+2*(NE+TableNum)
                                !输出低能区处理后的概率表
                                DO j=1,NE
                                    TotalPoint = LibIso4%FILE4_MT(n)%E(j)%TotalPoint
                                    IF(m(j).eq.2)THEN
                                        m(j)=3
                                        xss3(ptr4+1)=LibIso4%FILE4_MT(n)%Interpolation;xss3(ptr4+2)=m(j)
                                        xss3(ptr4+3)=-1;xss3(ptr4+4)=0;xss3(ptr4+5)=1
                                        xss3(ptr4+6)=0.5;xss3(ptr4+7)=0.5;xss3(ptr4+8)=0.5
                                        xss3(ptr4+9)=0;xss3(ptr4+10)=0.5;xss3(ptr4+11)=1
                                        ptr4=ptr4+2+2*m(j)
                                    ELSEIF(m(j).ne.2)THEN
                                        xss3(ptr4+1)=LibIso4%FILE4_MT(n)%Interpolation;xss3(ptr4+2)=m(j)
                                        DO i=1,m(j)
                                            xss3(ptr4+2+i)=angle2(i,j)
                                            xss3(ptr4+2+m(j)+i)=func(angle2(i,j),LibIso4%FILE4_MT(n)%E(j) &
                                                                         %Coeffi_a,TotalPoint)
                                        ENDDO
                                        ptr4=ptr4+2+2*m(j)
                                        pSumL(1)=0
                                        DO i=2,m(j)
                                            pSumL(i)=pSumL(i-1)+(func(angle2(i,j),LibIso4%FILE4_MT(n)%E(j)%Coeffi_a, &
                                                     TotalPoint)+func(angle2(i-1,j),LibIso4%FILE4_MT(n)%E(j)%Coeffi_a, &
                                                     TotalPoint))*(angle2(i,j)-angle2(i-1,j))/2
                                        ENDDO
                                        DO i=1,m(j)
                                            xss3(ptr4+i)=pSumL(i)/pSumL(m(j))
                                            xss3(ptr4-m(j)+i)=xss3(ptr4-m(j)+i)/pSumL(m(j))
                                        ENDDO         
                                        
                                    ENDIF
                                    ptr4=ptr4+m(j)
                                ENDDO
                                  !输出高能区概率表  
                                DO j=1,TableNum   
                                  DataNum=LibIso4%FILE4_MT(n)%Table(j)%DataNum
                                  xss3(ptr4+1)=LibIso4%FILE4_MT(n)%Table(j)%Interpolation
                                  xss3(ptr4+2)=DataNum
                                  DO l=1,DataNum
                                      xss3(ptr4+2+l)=LibIso4%FILE4_MT(n)%Table(j)%angular(l)
                                      xss3(ptr4+2+DataNum+l)=LibIso4%FILE4_MT(n)%Table(j)%probability(l)/LibIso4%FILE4_MT(n)%Table(j)%pSum(DataNum)
                                      xss3(ptr4+2+2*DataNum+l)=LibIso4%FILE4_MT(n)%Table(j)%pSum(l)/LibIso4%FILE4_MT(n)%Table(j)%pSum(DataNum)
                                  ENDDO
                                  ptr4=ptr4+2+3*DataNum
                                ENDDO
                                IF(LibIso4%FILE4_MT(n)%MT.eq.2)THEN
                                    ptr22=ptr4
                                ENDIF
                                IF(LibIso4%FILE4_MT(n)%MT.ne.2.and.ptr4-ptr3-total_number+1.gt.10.and.flag.eq.1) THEN
                                    flag=1
                                ENDIF
                                IF(n.lt.MT_total_number) THEN
                                    IF(LibIso4%FILE4_MT(n)%MT.eq.18.or.LibIso4%FILE4_MT(n)%MT.eq.2) THEN
                                        CONTINUE
                                    ELSE
                                        IF(flag.eq.1) THEN
                                            xss3(ptr44+jj)=ptr22-ptr3-total_number+1
                                            xss3(ptr44+jj+1)=ptr4-ptr3-total_number+1
                                            flag=0
                                        ELSE
                                            xss3(ptr44+jj+1)=ptr4-ptr3-total_number+1
                                        ENDIF
                                    ENDIF
                                ENDIF
                                
                                !DEALLOCATE(Loca)
                                 DEALLOCATE(angle2)
                                 DEALLOCATE(pSumL)
                                 DEALLOCATE(m)
                            ENDIF
                       ENDIF
                       EXIT
                   ENDIF
                   
               ENDDO
               count=ptr4-ptr3-total_number+1
           ENDDO
           ptr=ptr4
           mt18=0
           count=0
           print*,(mtr(i),i=1,mt_num)
           DO i=1,mt_num
               IF(mtr(i).lt.5.or.mtr(i).gt.91) THEN
                   count=count+1
               ELSE                   
                   DO j=1,MT_total_number5
                       DO k=1,MT_total_number6
                           IF(mtr(i).eq.LibIso5%FILE5_MT(j)%MT.and. &
                              mtr(i).eq.LibIso6%FILE6_MT(k)%MT) THEN
                               count=count+1
                           ENDIF
                       ENDDO
                   ENDDO
                ENDIF

           ENDDO
           
           jx=0
           ptr_ldlw=ptr   !LDLW
           jxs(10)=ptr_ldlw+1               !Table of energy distribution locators
           jxs(11)=ptr+mt_num-count+1       !Energy distributions
           ptr=ptr+mt_num-count    
           ptr_dlw=ptr
           DO i=1,mt_num
               m5=0
               m6=0
               IF(mtr(i).eq.18) THEN
                   !mt=18优先给file5处理
                   print*,MT_total_number5
                   DO j=1,MT_total_number5
                       print*,'LibIso5%FILE5_MT(j)%MT=',LibIso5%FILE5_MT(j)%MT
                       IF(LibIso5%FILE5_MT(j)%MT.eq.18) THEN
                           mt18=1
                           xss3(ptr_ldlw+i)=ptr-ptr_dlw+1
                           ALLOCATE(xss5(15000000))
                           !类似上面这种allocate有点浪费内存，但目前对计算速度影响不大，后期需要的时候可以再优化
                           CALL Self%File5DataProcess(FN2,j,xss5,ptr5,ptr,ptr_dlw)
                           PRINT*,'PTR5=',PTR5
                           DO k=1,ptr5-ptr
                               xss3(ptr+k)=xss5(ptr+k)
                           ENDDO
                           ptr=ptr5
                           DEALLOCATE(XSS5)
                           
                       ENDIF
                   ENDDO
                   !IF(mt18.eq.0) THEN
                   !    !mtt(i)=6
                   !ENDIF
               ELSE
                   IF(mtr(i).ge.51.and.mtr(i).le.91) THEN
                       DO j=1,MT_total_number5
                           IF(mtr(i).eq.LibIso5%FILE5_MT(j)%MT) THEN
                               m5=1
                           ENDIF   
                       ENDDO
                       
                       DO L=1,MT_total_number6
                            IF(mtr(i).eq.LibIso6%FILE6_MT(L)%MT)  THEN
                                m6=1
                            ENDIF
                       ENDDO
                       
                       
                       IF(m5.ne.1.and.m6.ne.1) THEN
                           DO k=1,100
                               IF(mtr(i).eq.LibIso%FILE3_MT(k)%MT) THEN
                                   xss3(ptr_ldlw+i)=ptr-ptr_dlw+1
                                   NP=LibIso%FILE3_MT(k)%NP
                                   AWR=LibIso4%FILE4_MT(1)%AWR
                                   xxx=(AWR+1)/AWR
                                   xss3(ptr+1)=0
                                   xss3(ptr+2)=3
                                   xss3(ptr+3)=ptr-ptr_dlw+10    !not determined
                                   xss3(ptr+4)=0
                                   xss3(ptr+5)=2
                                   xss3(ptr+6)=LibIso%FILE3_MT(k)%Energy_point(1)/mev
                                   xss3(ptr+7)=LibIso%FILE3_MT(k)%Energy_point(NP)/mev
                                   xss3(ptr+8)=1
                                   xss3(ptr+9)=1
                                   xss3(ptr+10)=xxx*abs(LibIso%FILE3_MT(k)%QI)/mev
                                   xss3(ptr+11)=1/(xxx*xxx)
                                   ptr=ptr+11
                                   EXIT
                                ENDIF
                           ENDDO
                       ELSEIF(m5.eq.1.and.m6.ne.1) THEN
                           DO j=1,MT_total_number5
                               IF(mtr(i).eq.LibIso5%FILE5_MT(j)%MT) THEN
                                   xss3(ptr_ldlw+i)=ptr-ptr_dlw+1
                                   ALLOCATE(xss5(15000000))
                                   CALL Self%File5DataProcess(FN2,j,xss5,ptr5,ptr,ptr_dlw)
                                   DO k=1,ptr5-ptr
                                       xss3(ptr+k)=xss5(ptr+k)
                                   ENDDO
                                   ptr=ptr5
                                   DEALLOCATE(xss5)
                               ENDIF
                           ENDDO
                       ELSEIF(m5.ne.1.and.m6.eq.1) THEN
                           DO j=1,MT_total_number6
                               IF(mtr(i).eq.LibIso6%FILE6_MT(j)%MT) THEN

                                   !mtt(i)=6
                                   !xss3(ptr_ldlw+i)=ptr-ptr_dlw+1
                                   ALLOCATE(xss6(20000000))
                                   CALL Self%File6DataProcess(FN2,j,xss6,ptr6,ptr,ptr_dlw,ptr_dlw2)
                                   xss3(ptr_ldlw+i)=ptr_dlw2-ptr_dlw+1

                                   DO k=1,ptr6-ptr
                                       xss3(ptr+k)=xss6(ptr+k)
                                   ENDDO
                                   ptr=ptr6
                                   DEALLOCATE(xss6)
                               ENDIF
                           ENDDO
                        ENDIF
                   ELSEIF (mtr(i).lt.51) THEN
                        DO j=1,MT_total_number5
                            
                            IF(mtr(i).eq.LibIso5%FILE5_MT(j)%MT) THEN
                                xss3(ptr_ldlw+i)=ptr-ptr_dlw+1
                                ALLOCATE(xss5(15000000))
                                CALL Self%File5DataProcess(FN2,j,xss5,ptr5,ptr,ptr_dlw)
                                DO k=1,ptr5-ptr
                                    xss3(ptr+k)=xss5(ptr+k)
                                ENDDO
                                ptr=ptr5
                                DEALLOCATE(xss5)
                            ENDIF
                        ENDDO
                        DO j=1,MT_total_number6
                            
                            IF(mtr(i).eq.LibIso6%FILE6_MT(j)%MT) THEN
                                !mtt(i)=6
                                !xss3(ptr_ldlw+i)=ptr-ptr_dlw+1
                                ALLOCATE(xss6(20000000))
                                CALL Self%File6DataProcess(FN2,j,xss6,ptr6,ptr,ptr_dlw,ptr_dlw2)
                                xss3(ptr_ldlw+i)=ptr_dlw2-ptr_dlw+1
                                DO k=1,ptr6-ptr
                                    xss3(ptr+k)=xss6(ptr+k)
                                ENDDO
                                ptr=ptr6
                                DEALLOCATE(xss6)
                            ENDIF
                        ENDDO
                   ENDIF
               ENDIF
           ENDDO
           
        !---------- 
        !unresolved
        !----------
        
        IF(LibIso2%MF2153%flag.eq.1) THEN    !如果存在不可分辨共振则处理
            jxs(23)=ptr+1
            xss3(ptr+1)=LibIso2%MF2153%num2
            xss3(ptr+2)=LibIso2%MF2153%num3  !不确定
            xss3(ptr+3)=2
            xss3(ptr+4)=LibIso2%MF2153%iinel
            xss3(ptr+5)=LibIso2%MF2153%iabso
            xss3(ptr+6)=LibIso2%MF2153%lssf
            ptr=ptr+6
            num2=LibIso2%MF2153%num2
            num3=LibIso2%MF2153%num3
            DO i=1,num2
                xss3(ptr+i)=LibIso2%MF2153%array(i,1)/mev
            ENDDO
            ptr=ptr+num2
            DO i=1,num2
                DO j=1,num3
                    IF(j.eq.1) THEN
                        xss3(ptr+j)=LibIso2%MF2153%array(i,j+1)
                    ELSEIF(j.gt.1) THEN
                        xss3(ptr+j)=LibIso2%MF2153%array(i,j+1)+xss3(ptr+j-1)
                    ENDIF
                ENDDO
                ptr=ptr+num3
                DO j=1,4*num3
                    xss3(ptr+j)=LibIso2%MF2153%array(i,1+num3+j)
                ENDDO
                ptr=ptr+4*num3
                IF(LibIso2%MF2153%lssf.eq.0) THEN
                    DO j=1,num3
                        xss3(ptr+j)=LibIso2%MF2153%array(i,1+5*num3+j)/mev
                    ENDDO
                ELSE
                    DO j=1,num3
                        xss3(ptr+j)=LibIso2%MF2153%array(i,1+5*num3+j)
                    ENDDO
                ENDIF
                ptr=ptr+num3
            ENDDO
            
        ENDIF
                
        !--------------    
        !delayed neutron   
        !---------------
        
        DO i=1,MT_total_number5
            count=0
            error=100
            tot=0
            IF(LibIso5%FILE5_MT(i)%MT.eq.455) THEN
                jxs(24)=ptr+1
                xss3(ptr+1)=LibIso1%FILE1%LNU455
                xss3(ptr+2)=0
                xss3(ptr+3)=LibIso1%FILE1%NP455
                DO j=1,LibIso1%FILE1%NP455
                    xss3(ptr+3+j)=LibIso1%FILE1%EintNuDelay(2*j-1)/mev
                    xss3(ptr+3+LibIso1%FILE1%NP455+j)=LibIso1%FILE1%EintNuDelay(2*j)
                ENDDO
                ptr=ptr+3+2*LibIso1%FILE1%NP455
                jxs(25)=ptr+1
                ALLOCATE(LibIso1%FILE1%DecayConstant(LibIso5%FILE5_MT(i)%Totle_NK),stat=status)
                DO j=1,LibIso5%FILE5_MT(i)%Totle_NK
                    xss3(ptr+1)=LibIso1%FILE1%DecayConstant(j)/shake
                    IF(LibIso5%FILE5_MT(i)%NK(j)%NR.eq.1.and.LibIso5%FILE5_MT(i)%NK(j)%TOTLE_PE.eq.2) THEN
                        xss3(ptr+2)=0
                    ELSE
                        PRINT*,'OTHER CONDITION_(MT=455)!'
                        PAUSE
                        STOP
                    ENDIF
                    xss3(ptr+3)=LibIso5%FILE5_MT(i)%NK(j)%NP
                    DO k=1,LibIso5%FILE5_MT(i)%NK(j)%NP
                        xss3(ptr+3+k)=LibIso5%FILE5_MT(i)%NK(j)%PE(k,1)/mev
                        xss3(ptr+3+LibIso5%FILE5_MT(i)%NK(j)%NP+k)=LibIso5%FILE5_MT(i)%NK(j)%PE(k,2)
                    ENDDO
                    ptr=ptr+3+2*LibIso5%FILE5_MT(i)%NK(j)%NP
                ENDDO
                ptr455=ptr
                jxs(26)=ptr455+1
                nxs(8)=LibIso5%FILE5_MT(i)%Totle_NK
                ptr=ptr+LibIso5%FILE5_MT(i)%Totle_NK   !提供455的位置指针空间
                jxs(27)=ptr+1
                DO j=1,LibIso5%FILE5_MT(i)%Totle_NK
                    xss3(ptr455+j)=ptr-ptr455-LibIso5%FILE5_MT(i)%Totle_NK+1
                    xss3(ptr+1)=0
                    xss3(ptr+2)=4
                    xss3(ptr+3)=ptr-ptr455-LibIso5%FILE5_MT(i)%Totle_NK+10
                    xss3(ptr+4)=0
                    xss3(ptr+5)=2
                    xss3(ptr+6)=LibIso5%FILE5_MT(i)%NK(1)%PE(1,1)/mev
                    xss3(ptr+7)=LibIso5%FILE5_MT(i)%NK(1)%PE(LibIso5%FILE5_MT(i)%NK(j)%NP,1)/mev
                    xss3(ptr+8)=1
                    xss3(ptr+9)=1
                    IF(LibIso5%FILE5_MT(i)%LF.eq.1) THEN
                        PRINT*,'455中出现LF=1!'
                        PAUSE
                        STOP
                    ELSEIF(LibIso5%FILE5_MT(i)%LF.eq.5) THEN
                        xss3(ptr+10)=0
                        xss3(ptr+11)=2
                        xss3(ptr+12)=LibIso5%FILE5_MT(i)%NK(1)%PE(1,1)/mev
                        xss3(ptr+13)=LibIso5%FILE5_MT(i)%NK(1)%PE(LibIso5%FILE5_MT(i)%NK(j)%NP,1)/mev
                        xss3(ptr+14)=16+count
                        error=LibIso5%FILE5_MT(i)%NK2(j)%g(2,1)                         
                        DO while(error.gt.40)      !从0到第二个能量点间距太大时，需扩展能量，tot为扩展的能量数
                            error=error*0.8409
                            tot=tot+1
                        ENDDO
                        num=LibIso5%FILE5_MT(i)%NK2(j)%TOTLE_E1int
                        
                        xss3(ptr+15)=16+count+2+3*(tot+num)
                        ptr=ptr+15
                        DO L=1,2   !输出两遍
                            
                            xss3(ptr+1)=LibIso5%FILE5_MT(i)%NK2(j)%intflag
                            xss3(ptr+2)=tot+num
                            ptr=ptr+2
                            xss3(ptr+1)=LibIso5%FILE5_MT(i)%NK2(j)%g(1,1)/mev
                            
                            DO k=1,tot
                                xss3(ptr+1+k)=(LibIso5%FILE5_MT(i)%NK2(j)%g(2,1)*0.8409**(tot-k+1))/mev                                
                            ENDDO
                            DO k=1,tot+1
                                xss3(ptr+tot+num+k)=(LibIso5%FILE5_MT(i)%NK2(j)%g(2,2)*0.1606*0.917**(tot+1-k))*mev
                            ENDDO
                            xss3(ptr+1+tot+num-1+1)=(LibIso5%FILE5_MT(i)%NK2(j)%g(2,2)*0.1606*0.917**(tot)*mev)
                            DO k=2,num
                                xss3(ptr+1+tot+k-1)=LibIso5%FILE5_MT(i)%NK2(j)%g(k,1)/mev
                                xss3(ptr+1+2*tot+num+k-1)=LibIso5%FILE5_MT(i)%NK2(j)%g(k,2)*mev
                            ENDDO
                            ptr=ptr+2*(tot+num)
                            xss3(ptr+1)=0
                            xss3(ptr+tot+num)=1
                            DO k=2,tot+num-1
                                xss3(ptr+k)=xss3(ptr+k-1)+(xss3(ptr-2*(tot+num)+k)-xss3(ptr-2*(tot+num)+k-1))*(xss3(ptr-(tot+num)+k)+xss3(ptr-(tot+num)+k))/2
                            ENDDO
                            ptr=ptr+tot+num
                        ENDDO
                        count=ptr-ptr455-LibIso5%FILE5_MT(i)%Totle_NK
                    ENDIF
                ENDDO
                         
            ENDIF
        ENDDO
    

! 处理File5 的数据
!CALL Self%File5DataProcess(FN)
!
!
! 处理File6 的数据
!CALL Self%File6DataProcess(FN)
!
! 概率表的输出处理
nxs(1)=ptr
jxs(22)=ptr
 WRITE(ENDFIO%ACE_OUT,'(8I9)')(nxs(i),i=1,16)
 WRITE(ENDFIO%ACE_OUT,'(8I9)')(jxs(i),i=1,32)
 WRITE(ENDFIO%ACE_OUT,'(4ES20.11)')(xss3(i),i=1,ptr)
 DEALLOCATE(XSS3)

END SUBROUTINE 

 
SUBROUTINE File5DataProcess_sub(Self,FN,iii,xss5,ptr5,from,ptr_dlw)
   !---------------------------------------------------------------
   ! 本程序主要处理File5 数据，生成ACE格式的概率
   !
   !---------------------------------------------------------------
  IMPLICIT NONE
  ! 全局变量
  CLASS(ACEPData_Bank):: Self
  CHARACTER(*) ,INTENT(IN):: FN 
  
  ! 局部变量
  TYPE(ENDFDataTYPE) :: LibIsoFile5
  
  INTEGER :: TOTLE_E1int,TOTLE_Eint,Totle_NK,Totle_NK2,MT,LF,LAW,NEP,NA,ND,LEP,LANG
  INTEGER,DIMENSION(:),ALLOCATABLE :: mtr
  INTEGER,DIMENSION(:),ALLOCATABLE :: mtt
  REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: xss5
  REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: temp     !临时存放数组(angle)
  REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: temp2    !临时存放数组(probility)
  REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: temp3    !临时存放数组(cummunitive probility)
  REAL(KIND=DBL)::fr(200,3)
  INTEGER :: MT_total_number,MT_total_number2,MT_total_number5,MT_total_number6
  INTEGER ::mt_num,status,max,NP,count,ix,mt18,ct
  INTEGER ::ptr5,ptr2,ptr3,from,ptr55,ptr_dlw
  INTEGER :: i,j,k,n,m,t,l,v,p,iii,xxx,y,x   !loop index
  REAL(KIND=DBL) :: delt,cx,cxx,fx,ex,et
  REAL(KIND=DBL):: error,x2(20),xm,angle(300),b,c,mev
  mt_num=0;x=1;y=1;max=65000;ptr5=from;mev=1000000
  ix=1;fx=0.8409;ex=40
  
  ! 读取FIle 5的数据
  CALL LibIsoFile5%ENDF_read(FN,5)
  
  ! 处理File 5的数据
!  ALLOCATE(xss5(500000))
                         
                         Totle_NK=LibIsoFile5%FILE5_MT(iii)%Totle_NK
                         DO k=1,Totle_NK
                                ptr5=ptr5+1
                                ct=ptr5
                                xss5(ptr5)=0
                                LF=LibIsoFile5%FILE5_MT(iii)%NK(k)%LF
                                xss5(ptr5+1)=LF              
                                IF(LF.eq.1.or.LF.eq.12) THEN
                                    xss5(ptr5+1)=4
                                ENDIF
                                
                                xss5(ptr5+2)=ptr5+10-ptr_dlw-1    !not determined
                                xss5(ptr5+3)=0
                                NP=LibIsoFile5%FILE5_MT(iii)%NP
                                xss5(ptr5+4)=NP
                                print*,'xss5(ptr5+4)=',xss5(ptr5+4)
                                DO l=1,NP 
                                   xss5(ptr5+4+l)=LibIsoFile5%FILE5_MT(iii)%NK(k)%PE(l,1)/mev
                                   xss5(ptr5+4+NP+l)=LibIsoFile5%FILE5_MT(iii)%NK(k)%PE(l,2)
                               ENDDO
                          
                                ptr5=ptr5+4+2*NP
                                IF(LF.eq.1) THEN
                                    !Totle_NK=LibIsoFile5%FILE5_MT(iii)%Totle_NK
                                    !DO k=1,Totle_NK
                                        TOTLE_Eint=LibIsoFile5%FILE5_MT(iii)%NK(k)%TOTLE_Eint  
                                        xss5(ptr5+1)=0
                                        xss5(ptr5+2)=TOTLE_Eint
                                        DO l=1,TOTLE_Eint
                                            xss5(ptr5+2+l)=LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint(l)/mev
                                            !xss5(ptr5+2+TOTLE_Eint+l)=ptr5-dl   !location not determined
                                        ENDDO
                                        ptr55=ptr5+2+TOTLE_Eint
                                        ptr5=ptr5+2+2*TOTLE_Eint
                                        DO l=1,TOTLE_Eint
                                            xss5(ptr55+l)=ptr5-ptr_dlw+1
                                            TOTLE_E1int=LibIsoFile5%FILE5_MT(iii)%NK(k)%TOTLE_E1int(l)
                                            count=TOTLE_E1int
                                            DO m=1,TOTLE_E1int-1
                                                delt=LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint_num(l)%g(m+1,1) &
                                                     -LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint_num(l)%g(m,1)
                                                IF(delt.gt.200000) THEN
                                                    count=count+4
                                                ENDIF
                                            ENDDO
                                            xss5(ptr5+1)=2
                                            xss5(ptr5+2)=count
                                            ALLOCATE(temp(count),stat=status)
                                            ALLOCATE(temp2(count),stat=status)
                                            ALLOCATE(temp3(count),stat=status)
                                            m=1
                                            ix=1
                                            DO WHILE(m.lt.count)
                                                temp(m)=LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint_num(l)%g(ix,1)/mev
                                                
                                                
                                                temp2(m)=LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint_num(l)%g(ix,2)*mev
                                                delt=LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint_num(l)%g(ix+1,1) &
                                                     -LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint_num(l)%g(ix,1)
                                                IF(delt.gt.200000) THEN
                                                    DO v=1,4
                                                        temp(m+v)=temp(m)+v*delt/5000000
                                                        temp2(m+v)=LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint_num(l)%g(ix,2)*exp((temp(m+v)- &
                                                                    LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint_num(l)%g(ix,1))*log(LibIsoFile5% &
                                                                     FILE5_MT(iii)%NK(k)%Eint_num(l)%g(ix+1,2)/LibIsoFile5%FILE5_MT(iii)% &
                                                                     NK(k)%Eint_num(l)%g(ix,2))/(LibIsoFile5%FILE5_MT(iii)%NK(k)% &
                                                                     Eint_num(l)%g(ix+1,1)-LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint_num(l)%g(ix,1)))
                                                        if(temp2(m+v).gt.1000000000000000) then    !少部分概率密度出现无穷大，可将其变成0,不影响结果
                                                            temp2(m+v)=0
                                                        endif
                                                    ENDDO
                                                    m=m+5
                                                    ix=ix+1
                                                    
                                                ELSE
                                                    m=m+1
                                                    ix=ix+1
                                                ENDIF                 
                                            ENDDO
                                            temp(count)=LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint_num(l)%g(TOTLE_E1int,1)/mev
                                            temp2(count)=LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint_num(l)%g(TOTLE_E1int,2)*mev
                                            temp3(1)=0
                                            !temp3(count)=1
                                            DO v=2,count
                                                temp3(v)=temp3(v-1)+(temp2(v-1)+temp2(v))*(temp(v)-temp(v-1))/2
                                            ENDDO
                                            DO v=1,count
                                                xss5(ptr5+2+v)=temp(v)
                                                xss5(ptr5+2+count+v)=temp2(v)/temp3(count)
                                                xss5(ptr5+2+2*count+v)=temp3(v)/temp3(count)
                                            ENDDO
                                            DEALLOCATE(temp)
                                            DEALLOCATE(temp2)
                                            DEALLOCATE(temp3)
                                            ptr5=ptr5+2+3*count  
                                        ENDDO
                                    !ENDDO
                                ELSE IF(LF.eq.5) THEN
                                    !Totle_NK=LibIsoFile5%FILE5_MT(iii)%Totle_NK
                                    !DO k=1,Totle_NK
                                        TOTLE_Eint=LibIsoFile5%FILE5_MT(iii)%NK(k)%TOTLE_Eint
                                        DO l=1,TOTLE_Eint
                                            xss5(ptr5+2+l)=LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint(l)/mev
                                            !xss5(ptr5+2+TOTLE_Eint+1)=404   !location not determined
                                        ENDDO
                                        ptr55=ptr5+2+TOTLE_Eint
                                        ptr5=ptr5+2+2*TOTLE_Eint
                                        DO l=1,TOTLE_Eint
                                            xss5(ptr55+l)=ptr5-ptr_dlw+2
                                            TOTLE_E1int=LibIsoFile5%FILE5_MT(iii)%NK(k)%Eint_num(l)%TOTLE_E1int
                                            xss5(ptr5+1)=0
                                            xss5(ptr5+2)=TOTLE_E1int
                                            PRINT*,'XSS5(I)=',xss5(ptr5+2)
                                            ALLOCATE(temp(TOTLE_E1int),stat=status)
                                            ALLOCATE(temp2(TOTLE_E1int),stat=status)
                                            ALLOCATE(temp3(TOTLE_E1int),stat=status)
                                            temp3(1)=0
                                            !temp3(TOTLE_E1int)=1
                                            DO m=2,TOTLE_E1int
                                                temp3(m)=temp3(m-1)+(temp2(m-1)+temp2(m))*(temp(m)-temp(m-1))/2
                                            ENDDO
                                            DO m=1,TOTLE_E1int
                                                xss5(ptr5+2+m)=temp(m)/mev
                                                xss5(ptr5+2+TOTLE_E1int+m)=temp2(m)*mev/temp3(TOTLE_E1int)
                                                xss5(ptr5+2+2*TOTLE_E1int+m)=temp3(m)/temp3(TOTLE_E1int)
                                            ENDDO
                                            ptr5=ptr5+2+3*TOTLE_E1int  
                                        ENDDO
                                    !ENDDO
                                ELSE IF(LF.eq.7.or.LF.eq.9) THEN
                                    CONTINUE
                                ELSE IF(LF.eq.12) THEN
                                    CONTINUE
                                ENDIF
                         ENDDO
                          !DEALLOCATE(xss5)
  WRITE(ENDFIO%OTAPE,'(4ES20.11)')(xss5(i),i=ct,ptr5)
  !WRITE(ENDFIO%OUTP,'(4ES21.11)')(xss5(i),i=1,ptr5)
  WRITE(ENDFIO%OUTP,*)'========='
END SUBROUTINE
   
   
    SUBROUTINE File6DataProcess_sub(Self,FN,xxx,xss6,ptr6,from,ptr_dlw,ptr_dlw2)
   !---------------------------------------------------------------
   ! 本程序主要处理File6数据，生成ACE格式的概率
   !
   !---------------------------------------------------------------
  IMPLICIT NONE
  ! 全局变量
  CLASS(ACEPData_Bank):: Self
  CHARACTER(*) ,INTENT(IN):: FN 
  
  ! 局部变量
  TYPE(ENDFDataTYPE) :: LibIsoFile6
  
  INTEGER :: TOTLE_E1int,TOTLE_Eint,Totle_NK,Totle_NK2,MT,LF,LAW,NEP2,NA,ND,LEP,LANG,Eint_num,TOTLE_YE
  INTEGER :: ZAP,IZAI,IVAR,NR
  INTEGER,DIMENSION(:),ALLOCATABLE :: mtr
  INTEGER,DIMENSION(:),ALLOCATABLE :: mtt
  !INTEGER,DIMENSION(:),ALLOCATABLE :: m2
  REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: xss6
  REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: temp     !临时存放数组(angle)
  REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: temp2    !临时存放数组(probility)
  REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: temp3    !临时存放数组(cummunitive probility)
  REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE :: coeff
  REAL(KIND=DBL)::fr(300,3)
  INTEGER :: MT_total_number,MT_total_number2,MT_total_number5,MT_total_number6
  INTEGER ::mt_num,status,max,NP,count,ix,mt18
  INTEGER ::ptr5,ptr2,ptr3,ptr6,ptr66,from,ptr_dlw,ptr_dlw2
  INTEGER :: i,j,k,n,m,t,l,v,p,iii,xxx,x,y,ii,nn,kk,jj   !loop index
  REAL(KIND=DBL) :: delt,cx,cxx,fx,ex,et,dy,f
  REAL(KIND=DBL):: error,x2(50),xm,angle(400),b,c,mev,mmm,y2(24),aco(7000),cprob(7000),cumm(7000)
  REAL(KIND=DBL):: probm,prob1,prob2,dm,ym,yt,xl,yl,check,dco,diff
  mt_num=0;x=1;y=1;max=65000;ptr6=from;mev=1000000
  ix=1;fx=0.8409;ex=40
  IZAI=1  !中子是1
  IVAR=0
  
  ! 读取File 6的数据
  CALL LibIsoFile6%ENDF_read(FN,6)
  
                          Totle_NK2 = LibIsoFile6%FILE6_MT(xxx)%Totle_NK     
                          DO k=1,Totle_NK2
                              TOTLE_Eint=LibIsoFile6%FILE6_MT(xxx)%NK(k)%TOTLE_Eint
                              LEP=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%LEP
                              TOTLE_YE=LibIsoFile6%FILE6_MT(xxx)%NK(k)%TOTLE_YE
                              LAW=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW
                              ZAP=LibIsoFile6%FILE6_MT(xxx)%NK(k)%ZAP
                              NR=LibIsoFile6%FILE6_MT(xxx)%NK(k)%NR
                              NP=LibIsoFile6%FILE6_MT(xxx)%NK(k)%NP
                              TOTLE_YE=LibIsoFile6%FILE6_MT(xxx)%NK(k)%TOTLE_YE
                              IF(ZAP.eq.0.or.ZAP-IZAI.gt.0.001_dbl) THEN  
                                  CONTINUE  !FILE6有多个nk，有些不处理
                              ELSE
                                  IF(TOTLE_YE.eq.1) THEN
                                      NP=NP-1
                                  ENDIF
                                  IF(NP.ge.2) THEN
                                      DO l=2,NP
                                          ii=LibIsoFile6%FILE6_MT(xxx)%NK(k)%YE(l,2)-LibIsoFile6%FILE6_MT(xxx)%NK(k)%YE(l-1,2)
                                          IF(ii.ne.0) THEN
                                              IVAR=1
                                          ENDIF
                                      ENDDO
                                  ENDIF
                                  ii=LibIsoFile6%FILE6_MT(xxx)%NK(k)%YE(1,2)-nint(LibIsoFile6%FILE6_MT(xxx)%NK(k)%YE(1,2))
                                  IF(IVAR.eq.0.and.abs(ii).gt.0.001_dbl) THEN
                                      IVAR=2
                                  ENDIF
                                  IF(IVAR.le.0) THEN
                                      CONTINUE
                                      !此时TYR正负以及数值需要变化，目前没搞明白，回头再改
                                  ELSE
                                     ptr6=ptr6+1
                                     xss6(ptr6)=0
                                     xss6(ptr6+1)=TOTLE_Eint
                                     ptr6=ptr6+1
                                     DO l=1,TOTLE_Eint
                                         xss6(ptr6+l)=LibIsoFile6%FILE6_MT(xxx)%NK(k)%YE(l,1)/mev
                                         xss6(ptr6+TOTLE_Eint+l)=LibIsoFile6%FILE6_MT(xxx)%NK(k)%YE(l,2)
                                     ENDDO
                                     ptr6=ptr6+2*TOTLE_Eint
                                  ENDIF     
                                  ptr_dlw2=ptr6
                                 IF(LAW.eq.1) THEN
                                    LANG=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%LANG
                                    Eint_num=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%Eint_num
                                    IF(LANG.eq.1) THEN
                                        !ptr6=ptr6+1
                                        !xss6(ptr6)=0
                                        !xss6(ptr6+1)=TOTLE_Eint
                                        !print*,TOTLE_Eint
                                        !pause
                                        !IF(TOTLE_Eint.gt.2) THEN
                                        !    DO l=1,TOTLE_Eint
                                        !        xss6(ptr6+1+l)=LibIsoFile6%FILE6_MT(xxx)%NK(k)%YE(l,1)/mev
                                        !        xss6(ptr6+1+TOTLE_Eint+l)=LibIsoFile6%FILE6_MT(xxx)%NK(k)%YE(l,2)
                                        !    ENDDO
                                        !    xss6(ptr6+1+2*TOTLE_Eint+1)=0
                                        !    xss6(ptr6+1+3*TOTLE_Eint)=1
                                        !    ptr6=ptr6+1+2*TOTLE_Eint
                                        !ENDIF
                                        xss6(ptr6+1)=0
                                        
                                        !IF(LAW.eq.1.and.LANG.ne.2) THEN
                                            xss6(ptr6+2)=61
                                   !     ELSEIF(LAW.eq.1.and.LANG.eq.2) THEN
                                        !    xss6(ptr6+2)=44
                                        !ENDIF
                                           
                                        xss6(ptr6+3)=ptr6-ptr_dlw+14   
                                        xss6(ptr6+4)=TOTLE_YE
                                        xss6(ptr6+5)=TOTLE_Eint
                                        xss6(ptr6+6)=TOTLE_YE
                                        xss6(ptr6+7)=TOTLE_Eint
                                        ptr6=ptr6+7
                                        
                                      
                                        DO l=1,TOTLE_Eint 
                                            xss6(ptr6+l)=LibIsoFile6%FILE6_MT(xxx)%NK(k)%YE(l,1)/mev
                                            xss6(ptr6+TOTLE_Eint+l)=1
                                        ENDDO
                                        ptr6=ptr6+2*TOTLE_Eint
                                        xss6(ptr6+1)=0
                                        xss6(ptr6+2)=Eint_num
                                        DO l=1,Eint_num
                                            xss6(ptr6+2+l)=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%Eint(l)/mev
                                        ENDDO
                                        ptr6=ptr6+2+Eint_num
                                        ptr3=ptr6
                                        ptr6=ptr6+Eint_num   !给每一个入射能量的位置指针提供空间
                                        DO L=1,Eint_num
                                            xss6(ptr3+L)=ptr6-ptr_dlw+1
                                            ND=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%ND
                                            NEP2=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%NEP(l)
                                            NA=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%NA(l)
                                            
                                            xss6(ptr6+1)=LEP+10*ND
                                            xss6(ptr6+2)=NEP2
                                            DO m=1,NEP2
                                                xss6(ptr6+2+m)=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%E_num(l)%b(m,1)/mev
                                                xss6(ptr6+2+NEP2+m)=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%E_num(l)%b(m,2)*mev
                                            ENDDO
                                            IF(NEP2.eq.2) THEN
                                                xss6(ptr6+2+2*NEP2+1)=0
                                                xss6(ptr6+2+3*NEP2)=1
                                            ELSEIF(NEP2.gt.2) THEN
                                               xss6(ptr6+2+2*NEP2+1)=0
                                               xss6(ptr6+2+3*NEP2)=1
                                               DO m=2,NEP2-1
                                                   xss6(ptr6+2+2*NEP2+m)=xss6(ptr6+2+2*NEP2+m-1)+(xss6(ptr6+2+NEP2+m)+xss6(ptr6+2+NEP2+m-1))*(xss6(ptr6+2+m)-xss6(ptr6+2+m-1))/2
                                               ENDDO
                                            ENDIF
                                            ptr6=ptr6+2+3*NEP2
                                            ptr2=ptr6
                                            ptr6=ptr6+NEP2    !给入射能量的每一个出射能量的位置指针提供空间
                                            
                                            DO m=1,NEP2
                                                ALLOCATE(coeff(NA))
                                                DO v=1,NA
                                                    coeff(v)=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%E_num(l)%b(m,2+v)/LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%E_num(l)%b(m,2)
                                                ENDDO
                                                xss6(ptr2+m)=ptr6-ptr_dlw+1  
                                                xss6(ptr6+1)=2
                                                !ptr66=ptr6+2
                                                !ptr6=ptr6+2
                                                ii=0
                                               ! DO p=1,NA
                                                    count=1
                                                    cumm(1)=0
                                                   x2(1)=1;x2(2)=-1;c=x2(1);b=x2(2);t=2
                                                   y2(1)=func(x2(1),coeff,NA)
                                                   y2(2)=func(x2(2),coeff,NA)
                                                   DO WHILE (t.gt.0)
                                                      dy=0
                                                      IF (t.gt.1.and.t.lt.24) THEN
                                                         dm=x2(t-1)-x2(t)
                                                         xm=(x2(t-1)+x2(t))/2
                                                         
                                                         IF (xm.gt.x2(t).and.xm.lt.x2(t-1)) THEN
                                                            ym=(y2(t-1)+y2(t))/2
                                                            yt=func(xm,coeff,NA)
                                                            error=0.0002_DBL*abs(yt)+1.0_DBL/1000000
                                                            dy=abs(yt-ym)
                                                            IF (dm.gt.10.0_dbl) dy=2*error
                                                            IF (ym.ne.0.and.yt/ym.gt.1.0_dbl*5) dy=2*error
                                                            IF (ym.ne.0.and.yt/ym.lt.1.0_dbl/5) dy=2*error
                                                         ENDIF
                                                      ENDIF
                                                      IF (dy.gt.error) THEN
                                                         ! not converged.
                                                         ! add midpoint to stack and continue.
                                                         t=t+1
                                                         x2(t)=x2(t-1)
                                                         y2(t)=y2(t-1)
                                                         x2(t-1)=xm
                                                         y2(t-1)=yt
                                                      ELSE
                                                         ! converged.
                                                         ! use top point in stack.
                                                         ii=ii+1
                                                         aco(ii)=x2(t)
                                 
                                                         cprob(ii)=y2(t)
                                                         IF (ii.gt.1) cumm(ii)=cumm(ii-1)+(x2(t)-xl)*(y2(t)+yl)/2
                                                         xl=x2(t)
                                                         yl=y2(t)
                                                         t=t-1
                                                      ENDIF
                                                   ENDDO
                                                   nn=ii-1
                                                   !thin  
                                                   i=1
                                                   ii=1
                                                   aco(ii)=-1
                                                   DO WHILE (i.lt.nn-1)
                                                      check=0
                                                      dco=0
                                                      j=i+1
                                                      DO WHILE (j.lt.nn+1.and.check.le.0.and.dco.le.1.0_dbl)
                                                         j=j+1
                                                         jj=j-1
                                                         dco=aco(j)-aco(i)
                                                         IF (dco.le.1.0_dbl) THEN
                                                            kk=i
                                                            DO WHILE (k.lt.j-1.and.check.le.0)
                                                               kk=kk+1
                                                               f=(aco(j)-aco(kk))/dco
                                                               error=f*cprob(i)+(1-f)*cprob(j)
                                                               diff=1.0_dbl/10000000+0.002_dbl*cprob(kk)
                                                               check=abs(error-cprob(kk))-diff
                                                            ENDDO
                                                         ENDIF
                                                      ENDDO
                                                      IF (check.gt.0.or.dco.gt.1.0_dbl) THEN
                                                         i=jj
                                                         ii=ii+1
                                                         aco(ii)=aco(i)
                                                         cprob(ii)=cprob(i)
                                                         IF (ii.gt.1) cumm(ii)=cumm(ii-1)&
                                                           +(aco(ii)-aco(ii-1))*(cprob(ii)+cprob(ii-1))/2
                                                      ENDIF
                                                   ENDDO
                                                   i=nn+1
                                                   ii=ii+1
                                                   aco(ii)=aco(i)
                                                   cprob(ii)=cprob(i)
                                                   cumm(ii)=cumm(ii-1)+(aco(ii)-aco(ii-1))*(cprob(ii)+cprob(ii-1))/2
                                 
                                                   IF(ii.eq.2) THEN
                                                       ii=3
                                                       xss6(ptr6+2)=ii
                                                       ptr6=ptr6+2
                                                       !xss6(ptr66)=ptr6+1
                                                       DO v=1,3
                                                           xss6(ptr6+v)=-2.0_dbl+v*1.0_dbl
                                                           xss6(ptr6+3+v)=0.5_dbl
                                                           xss6(ptr6+6+v)=-0.5_dbl+v*0.5_dbl
                                                       ENDDO
                                                   ELSE IF(ii.gt.2) THEN
                                                       xss6(ptr6+2)=ii
                                                       ptr6=ptr6+2
                                                       !xss6(ptr66)=ptr6+1
                                                       DO v=1,ii
                                                           xss6(ptr6+v)=aco(v)
                                                           xss6(ptr6+ii+v)=cprob(v)/cumm(ii)
                                                       ENDDO
                                                       xss6(ptr6+2*ii+1)=0
                                                       xss6(ptr6+3*ii)=1
                                                       DO v=2,ii-1
                                                           xss6(ptr6+2*ii+v)=xss6(ptr6+2*ii+v-1)+(cprob(v)/cumm(ii)+cprob(v-1)/cumm(ii))*(aco(v)-aco(v-1))/2
                                                       ENDDO
                                                   ENDIF
                                               ! ENDDO
                                               
                                               
                                                ptr6=ptr6+3*ii
                                                DEALLOCATE(coeff)
                                                
                                            ENDDO
                                        ENDDO
                                        
                                    ELSEIF(LANG.eq.2) THEN
                                        ptr6=ptr6+1
                                        xss6(ptr6)=0
                                        xss6(ptr6+1)=44
                                        xss6(ptr6+2)=ptr6-ptr_dlw+5+2*TOTLE_Eint
                                        xss6(ptr6+3)=0
                                        xss6(ptr6+4)=TOTLE_Eint
                                        DO l=1,TOTLE_Eint
                                            xss6(ptr6+4+l)=LibIsoFile6%FILE6_MT(xxx)%NK(k)%YE(l,1)/mev
                                            xss6(ptr6+4+TOTLE_Eint+l)=1  !有时候不全是1，回头再确定具体数值
                                        ENDDO
                                        ptr6=ptr6+4+2*TOTLE_Eint
                                        xss6(ptr6+1)=0
                                        xss6(ptr6+2)=Eint_num
                                        ptr6=ptr6+2
                                        DO l=1,Eint_num
                                            xss6(ptr6+l)=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%Eint(l)/mev
                                        ENDDO
                                        ptr6=ptr6+Eint_num
                                        ptr2=ptr6     !给每个入射能量提供位置空间
                                        ptr6=ptr6+Eint_num
                                        
                                        DO l=1,LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%Eint_num
                                            xss6(ptr2+l)=ptr6-ptr_dlw+1
                                            !xss6(ptr3+l)=ptr6
                                            ND=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%ND
                                            NEP2=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%NEP(l)
                                            xss6(ptr6+1)=LEP+10*ND
                                            fr(1:NEP2,:)=LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%E_num(l)%b(1:NEP2,:)
                                            cx=fr(2,1)*fr(1,2)
                                            t=2
                                            DO WHILE(NEP2.gt.2)
                                                cxx=cx+(fr(t+1,1)-fr(t,1))*fr(t,2)
                                                IF(abs(cxx/fr(t+1,1)**1.5-cx/fr(t,1)**1.5)>cx/fr(t,1)**1.5/50)exit                                            
                                                fr(1,2)=cxx/fr(t+1,1)                                          
                                                DO m=3,NEP2                                        
                                                    fr(m-1,:)=fr(m,:)
                                                ENDDO
                                                cx=cxx
                                                NEP2=NEP2-1
                                                t=t+1
                                            ENDDO
                                            et=fr(2,1)
                                            DO WHILE(et.gt.ex)
                                                et=et*fx
                                                DO t=NEP2,2,-1
                                                    fr(t+1,:)=fr(t,:)
                                                ENDDO
                                                fr(2,1)=et
                                                fr(2,2)=(1-fx*sqrt(fx))*fr(1,2)/(1-fx)
                                                fr(2,3)=fr(1,3)
                                                fr(1,2)=fr(1,2)*sqrt(fx)
                                                NEP2=NEP2+1
                                            ENDDO
                                            xss6(ptr6+2)=NEP2
                                            ptr6=ptr6+2
                                            DO m=1,NEP2
                                                xss6(ptr6+m)=fr(m,1)/mev                                                       ! E_OUT
                                                xss6(ptr6+NEP2+m)=fr(m,2)*mev                                                  ! Probability density function
                                                xss6(ptr6+3*NEP2+m)=fr(m,3)       ! Precompound fraction r
                                                xss6(ptr6+4*NEP2+m)=kalbach_a(1,1,LibIsoFile6%FILE6_MT(xxx)%IZAI, &            ! Angular distribution slope value a
                                                                         LibIsoFile6%FILE6_MT(xxx)%NK(k)%LAW1%Eint(l),fr(m,1)) 
                                               ! print*,xss6(ptr6+NEP2+m)
                                            ENDDO
                                            !pause
                                            xss6(ptr6+2*NEP2+1)=0
                                            !xss6(ptr6+3*NEP2)=1
                                            IF(LEP.eq.1) THEN                                           
                                                DO m=2,NEP2
                                                    xss6(ptr6+2*NEP2+m)=xss6(ptr6+2*NEP2+m-1)+fr(m,2)*(fr(m,1)-fr(m-1,1))              ! Cumulative density function
                                                ENDDO
                                            ELSEIF(LEP.eq.2) THEN
                                                DO m=2,NEP2
                                                    xss6(ptr6+2*NEP2+m)=xss6(ptr6+2*NEP2+m-1)+(fr(m,2)+fr(m-1,2))*(fr(m,1)-fr(m-1,1))/2
                                                ENDDO
                                            ENDIF
                                            DO m=1,NEP2                                             
                                               xss6(ptr6+2*NEP2+m)=xss6(ptr6+2*NEP2+m)/xss6(ptr6+3*NEP2)      !归一化
                                            ENDDO
                                            ptr6=ptr6+5*NEP2
                                        ENDDO
                                    ENDIF
                                 ELSEIF(LAW.eq.2) THEN
                                     ptr6=ptr6+1
                                     xss6(ptr6)=0
                                     xss6(ptr6+1)=66
                                     print*,'CAN NOT PROCESS LAW=2'
                                     PAUSE
                                     STOP
                                 ELSEIF(LAW.eq.5) THEN
                                     ptr6=ptr6+1
                                     xss6(ptr6)=0
                                     xss6(ptr6+1)=66
                                     print*,'CAN NOT PROCESS LAW=5'
                                     PAUSE
                                     STOP
                                 ELSEIF(LAW.eq.6) THEN
                                     ptr6=ptr6+1
                                     xss6(ptr6)=0
                                     xss6(ptr6+1)=66
                                    print*,'CAN NOT PROCESS LAW=6'
                                    PAUSE
                                     STOP
                                 ELSEIF(LAW.eq.7) THEN
                                     ptr6=ptr6+1
                                     xss6(ptr6)=0
                                     xss6(ptr6+1)=67
                                    print*,'CAN NOT PROCESS LAW=7'
                                    PAUSE
                                     STOP
                                 ENDIF
                              ENDIF
                          ENDDO
  
  
  END SUBROUTINE
   
   
  
  !
  !
  SUBROUTINE Legendre(x,legen,n)    !勒让德递推式
      IMPLICIT NONE
      INTEGER::n,k,i
      REAL(KIND=DBL)::x,g,h
      REAL(KIND=DBL)::legen(*)
      legen(1)=1
      legen(2)=x
      k=n-1
      DO i=1,k
          g=x*legen(i+1)
          h=g-legen(i)
          legen(i+2)=h+g-h/(i+1)
      ENDDO
      RETURN
  END SUBROUTINE Legendre
  
  
  REAL FUNCTION func(x,a,n)         !x是角度，a是勒让德系数，n是阶数，func算出概率
      IMPLICIT NONE
      REAL(KIND=DBL),INTENT(IN):: x
      INTEGER ::i
      INTEGER,INTENT(IN)::n
      REAL(KIND=DBL):: a(n)
      REAL(KIND=DBL)::p(100)
      func=0.5
      call Legendre(x,p,n)
      DO i=1,n
          func=func+(2*i+1)*a(i)*p(i+1)/2
      ENDDO
  END FUNCTION func
  
  
   REAL(DBL) function kalbach_a(aaa,bbb,zzaa,e,e_) 
       !求kalbach-mann分布的a
       !a+A->C->b+B
       !输入a,b的质量数，A的质子数与质量数，a的入射能量和b在质心系下的出射能量
       INTEGER::aaa,bbb,zzaa
       REAL(DBL)::e,e_
      REAL(DBL)::aa,ab,ac,na,nb,nc,za,zb,zc,sa,sb
      REAL(DBL)::ecm,ea,eb,x1,x3,fa,fb,a
      REAL(DBL),PARAMETER::third=1/3d0
      REAL(DBL),PARAMETER::twoth=2/3d0
      REAL(DBL),PARAMETER::fourth=4/3d0
      REAL(DBL),PARAMETER::c1=15.68d0
      REAL(DBL),PARAMETER::c2=-28.07d0
      REAL(DBL),PARAMETER::c3=-18.56d0
      REAL(DBL),PARAMETER::c4=33.22d0
      REAL(DBL),PARAMETER::c5=-0.717d0
      REAL(DBL),PARAMETER::c6=1.211d0
      REAL(DBL),PARAMETER::s2=2.22d0
      REAL(DBL),PARAMETER::s3=8.48d0
      REAL(DBL),PARAMETER::s4=7.72d0
      REAL(DBL),PARAMETER::s5=28.3d0
      REAL(DBL),PARAMETER::et1=130d0
      REAL(DBL),PARAMETER::et3=41d0
      REAL(DBL),PARAMETER::half=0.5d0
      REAL(DBL),PARAMETER::b1=0.04d0
      REAL(DBL),PARAMETER::b2=1.8d-6
      REAL(DBL),PARAMETER::b3=6.7d-7
      REAL(DBL),PARAMETER::d1=9.3d0
      REAL(DBL),PARAMETER::tomev=1.d-6  
       IF (aaa.eq.0) aaa=1         
       IF (zzaa.eq.12000) zzaa=12024
       aa=mod(zzaa,1000)
       za=int(zzaa/1000)
       ac=aa+mod(aaa,1000)
       zc=za+int(aaa/1000)
       ab=ac-mod(bbb,1000)
       zb=zc-int(bbb/1000)
       na=nint(aa-za)
       nb=nint(ab-zb)
       nc=nint(ac-zc)
       sa=c1*(ac-aa)+c2*((nc-zc)*(nc-zc)/ac-(na-za)*(na-za)/aa)+c3*(ac**twoth-aa**twoth)&
       +c4*((nc-zc)*(nc-zc)/ac**fourth-(nb-zb)*(nb-zb)/ab**fourth)&
       +c5*(zc*zc/ac**third-za*za/aa**third)+c6*(zc*zc/ac-za*za/aa)
       sb=c1*(ac-ab)+c2*((nc-zc)*(nc-zc)/ac-(nb-zb)*(nb-zb)/ab)+c3*(ac**twoth-ab**twoth)&
       +c4*((nc-zc)*(nc-zc)/ac**fourth-(nb-zb)*(nb-zb)/ab**fourth)&
       +c5*(zc*zc/ac**third-zb*zb/ab**third)+c6*(zc*zc/ac-zb*zb/ab)
       ecm=aa*e/ac
      ea=ecm*tomev+sa
      eb=tomev*e_*ac/ab+sb
      x1=eb
      IF (ea.gt.et1) x1=et1*eb/ea
      x3=eb
      IF (ea.gt.et3) x3=et3*eb/ea
      fa=1
      fb=1
      IF (bbb.eq.1) fb=half
      a=b1*x1+b2*x1**3+b3*fa*fb*x3**4
      kalbach_a=a
       return
       END function
  




END MODULE  ACEPModule
    