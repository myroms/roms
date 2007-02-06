!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name: MCT_1_0_12 $ 
!-----------------------------------------------------------------------
subroutine mpe_logging_init(myid,a,b,c,d,e,f,g,h,i,j,k,l,m,n)

implicit none
include 'mpif.h'

integer myid,a,b,c,d,e,f,g,h,i,j,k,l,ier
integer m,n,p,q,r,s,t,u,v
integer m1,m2,m3,m4,m5,m6
      
integer MPE_Log_get_event_number
integer MPE_Describe_state

#ifdef MPE
call MPE_init_log

a = MPE_Log_get_event_number()
b = MPE_Log_get_event_number()
c = MPE_Log_get_event_number()
d = MPE_Log_get_event_number()
e = MPE_Log_get_event_number()
f = MPE_Log_get_event_number()
g = MPE_Log_get_event_number()
h = MPE_Log_get_event_number()
i = MPE_Log_get_event_number()
j = MPE_Log_get_event_number()
k = MPE_Log_get_event_number()
l = MPE_Log_get_event_number()
m = MPE_Log_get_event_number()
n = MPE_Log_get_event_number()

if (myid .eq. 0) then

!        Description of init
ier = MPE_Describe_state( a, b, "MCT_Init", "coral")

         ier = MPE_Describe_state( c, d, "GSMap_ini", "maroon" )

         ier = MPE_Describe_state( e, f, "aV_ini", "aquamarine" )

         ier = MPE_Describe_state( g, h, "Rout_ini", "magenta" )

         ier = MPE_Describe_state( i, j, "MCT_Send", "yellow"      )

         ier = MPE_Describe_state( k, l, "MCT_Recv", "green"       )

         ier = MPE_Describe_state( m, n, "cleanup", "navy"       )

!        ier = MPE_Describe_state( p, q, "other", "firebrick"       )
!        ier = MPE_Describe_state( m1, m2, "other", "lavender"       )
!        ier = MPE_Describe_state( m3, m4, "other", "cyan"       )
!        ier = MPE_Describe_state( m5, m6, "other", "gray"       )

      endif
#endif
      return
      end
