MODULE MOD_ANGLE
	REAL*8,SAVE :: pi,rd
END MODULE

!*********************

MODULE MOD_staticvar
	REAL*8,SAVE :: cts,cto,ctscto
	REAL*8,SAVE :: ddb,ddf,dob,dof,dso
	REAL*8,SAVE :: ko,ks,sdb,sdf
	REAL*8,SAVE :: sob,sof,sumint
	REAL*8,ALLOCATABLE,SAVE :: sb(:),sf(:),vb(:),vf(:),w(:)
	REAL*8,ALLOCATABLE,SAVE :: m(:),m2(:),att(:),sigb(:),rinf(:),lidf(:)
END MODULE

!*********************

MODULE MOD_output_PROSPECT
	REAL*8,ALLOCATABLE,SAVE :: LRT(:,:),rho(:),tau(:)
END	MODULE

!*********************

MODULE MOD_flag_util
	LOGICAL flag(7),init_completed,init_completed0
	LOGICAL delta_geom,delta_lai,delta_hot,delta_leaf,delta_skyl,delta_soil,delta_lidf
	REAL*8,ALLOCATABLE,SAVE :: lidf_old(:),rsoil_old(:)
	REAL*8,SAVE :: N_old,Cab_old,Car_old,Cbrown_old,Cw_old,Cm_old
	REAL*8,SAVE :: LIDFa_old,LIDFb_old,lai_old,q_old,skyl_old
	REAL*8,SAVE :: tts_old,tto_old,psi_old
END	MODULE

!*********************

MODULE MOD_SAIL
	REAL*8,SAVE :: tss,too,tsstoo
	REAL*8,ALLOCATABLE,SAVE :: tso(:),tsd(:),tdd(:),tdo(:),rsd(:),rdd(:),rso(:),rdo(:)
	REAL*8,ALLOCATABLE,SAVE :: rsos(:),rsod(:),rddt(:),rsdt(:),rdot(:),rsodt(:),rsost(:),rsot(:)
END MODULE

!*********************