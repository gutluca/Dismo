C=======================================================================
C  Disease Impact and Severity Module  --  DISMO
C  Gustavo de Angelo Luca, Izael Martins Fattori Jr, Willingthon Pavan, Fabio Ricardo Marin
C  Luiz de Queiroz College of Agriculture (ESALQ),
C  University of Sao Paulo, Piracicaba, Brazil
C
C-----------------------------------------------------------------------
C  REVISION HISTORY
C  07/17/2023 Written.
C  11/15/2024 Revised.
C  05/20/2025 Fungicide logic.
C  08/10/2025 Logic/robustness fixes (cohorts, IR, LAF, LWD, cum. LAI)
C  11/05/2025 Write func for "DISMO.OUT" improved
C  11/17/2025 Virtual lesions factor added.
C  12/09/2025 Severity calculation in output file
C  03/05/2026 Parameter file lookup in the working directory.
C  03/10/2026 Moved DISMO.for from Plant\CROPGRO to Plant\Generic-Pest
C  06/25/2026 Improved output file formatting
C  07/17/2026 Added defoliation/senescence logic
C  07/28/2026 Added monocyclic disease support (NCYCLE parameter: M/P)
C  09/02/2026 Added pre-plant environmental inoculum reconstruction, primary inoculum
C             build-up and disease start based on weather favourability (FAV_THR) or observed onset (DAE_MIN).
C-----------------------------------------------------------------------
      SUBROUTINE DISEASE_LEAF (DYNAMIC,
     &    CONTROL, ISWITCH, Tmin, Tmax, RH, LAI_TOTAL,    ! Input
     &    WTLF, SLDOT,
     &    YRDOY, YREMRG, NVEG0, YREND,                    ! Input
     &    DISEASE_LAI, VIRTUAL_PHOTO_FACTOR,              ! Output
     &    DISEASE_SEN_RATE)                               ! Output
C-----------------------------------------------------------------------
      USE ModuleDefs
      IMPLICIT NONE
      EXTERNAL F_IR, F_DS, F_CANSPO, F_LR, F_LS, F_LAF, F_LEXP,
     &         F_LAR, F_PPSR_POP, APPLY_FUNGICIDE, CALC_DVIP,
     &         READ_DISEASE_PARAMETERS, F_VIRTUAL_LESIONS,
     &         F_SEVERITY, GETLUN, F_DEFOLIATION,
     &         DISMO_PRESEASON, DISMO_UPDATE_ENVIRONMENT,
     &         DISMO_SEASON_RESET, DISMO_CLEAR_SLOT, WARNING
      
      SAVE
C-----------------------------------------------------------------------
      INTEGER, PARAMETER :: MAXDAYS      = 250
      INTEGER, PARAMETER :: MAXPRESEASON = 400
      REAL,    PARAMETER :: EPS          = 1.0E-6
 
C  Named columns of the cohort array 
      INTEGER, PARAMETER :: C_LES  = 1   ! lesions m-2 created that day
      INTEGER, PARAMETER :: C_LATP = 2   ! latency progress (0..1)
      INTEGER, PARAMETER :: C_INFF = 3   ! infectious flag (0/1)
      INTEGER, PARAMETER :: C_AGE  = 4   ! relative lesion age (0..1)
 
C----- Dummy arguments -------------------------------------------------
      INTEGER DYNAMIC
      REAL    Tmin, Tmax, RH, LAI_TOTAL
      REAL    WTLF, SLDOT
      INTEGER YRDOY, YREMRG, NVEG0, YREND
      REAL    DISEASE_LAI, VIRTUAL_PHOTO_FACTOR, DISEASE_SEN_RATE
 
      TYPE (ControlType) CONTROL
      TYPE (SwitchType)  ISWITCH
 
C----- Calibrated parameters
      REAL    LESION_S, KVERHULST, RVERHULST
      REAL    YMAX, COF_A, COF_B
      REAL    TMIN_G, TOT_G, TMAX_G
      REAL    TMIN_D, TOT_D, TMAX_D
      REAL    LDMIN, LESIONAGEOPT, LESLIFEMAX
      REAL    BETA, RRDS
      REAL    SRC_HALF, FAV_THR, NDS, DEP_FRAC
      INTEGER DAE_MIN
      LOGICAL DAE_MIN_PRESENT, IS_MONOCYCLIC
      CHARACTER(LEN=1) NCYCLE
 
C----- Fixed model constants.
      REAL    LAI_MIN_START
      REAL    SPOR_DECAY, SRC_SURV, SEC_DECAY, FUNG_EFFICIENCY
      INTEGER FUNG_RES_D, FUNG_BUF_D, DVIP_THR, INOC_LOOKBACK
      LOGICAL USE_FUNGICIDE, USE_WTH_RH

C----- Environmental state -----------
      REAL    SOURCE_PRESSURE
      REAL    SEC_SPORE_CLOUD, SEC_SPORES_PENDING
      REAL    FAV_SUM, PRI_POOL
      LOGICAL PRI_RELEASED
      LOGICAL PRESEASON_DONE
      INTEGER PRESEASON_COUNT
      INTEGER PRESEASON_DATE(MAXPRESEASON)
      REAL    PRESEASON_RH(MAXPRESEASON)
      REAL    PRESEASON_LWD(MAXPRESEASON)
      REAL    PRESEASON_FT(MAXPRESEASON)
      REAL    PRESEASON_FAV(MAXPRESEASON)
 
C----- Epidemic state (season scope) -----------------------------------
C  PEND_AREA : leaf area claimed by live cohorts that has not yet turned
C              necrotic. Carried one day, like PREV_IS.
      REAL    PEND_AREA
      REAL    LAI_PEAK_SEASON, CUM_NECROTIC, PREV_IS
      REAL    SEVERITY_PCT
      REAL    ESP_LAT_HIST(MAXDAYS,4)
      REAL    ADMITTED_AREA(MAXDAYS)
      INTEGER DAE, PLANT_LIVE, N_ACTIVE_COH
      INTEGER SLOT_DAE(MAXDAYS)
 
C----- Fungicide / risk index state (season scope) ---------------------
      INTEGER DVIP_pts(7), idx, SUM7
      INTEGER BufferDays, ResidualDays, NSprays
      LOGICAL FungActive
 
C----- Output state ----------------------------------------------------
      INTEGER LUN_OUT
      LOGICAL HDR_DONE, LUN_OPEN, OVERFLOW_WARNED
 
C----- Daily working variables --------------------------------
      REAL    T, LWD, FT, FT_D, FT_G, DAILY_IP
      REAL    RH_OUT, LWD_OUT, FT_OUT, FAV_OUT
      REAL    IR, FSS, LR, LA, LAF, Lesion_Rate
      REAL    EXPN, EXPN_PREV, AREA_GROW, PEND_TODAY, CLAIM_TODAY
      REAL    HEALTH_LAI, LAI_SUSC, LAI_AVAIL, NEW_LOSS_TODAY
      REAL    IS, ESP_INOC_SEC, INF_AREA_K, LES_COUNT_PREV
      REAL    POT_SPO_PER_AREA, PS_K
      REAL    F_SOURCE
      REAL    PRI_CLOUD, SEC_CLOUD, CLOUD_TOTAL
      REAL    DS_TOTAL, DS_PRI, DS_SEC, LS_TODAY
      REAL    SEVFRAC, FVL
      INTEGER k, DAE_IDX, DVIP_today
      INTEGER IYEAR, IDOY, IDAP, I_PRE, IOS
      LOGICAL DEPOSITION_OK, EPI_ACTIVE
      CHARACTER(LEN=78) MSG(4)
 
C***********************************************************************
C  RUNINIT -- once per run
C***********************************************************************
      IF (DYNAMIC .EQ. RUNINIT) THEN
 
          USE_WTH_RH = .TRUE.
 
C-----------------------------------------------------------------------
C  Fixed model constants. 
C-----------------------------------------------------------------------
          SPOR_DECAY = 0.7937
          SEC_DECAY  = 0.7937
          SRC_SURV   = 0.98

C  Minimum LAI for deposition.
          LAI_MIN_START = 0.5
 
C  Pre-plant weather replay window (d) used to build initial inoculum.
          INOC_LOOKBACK = 60
 
C  Fungicide block -- move to the FILEX management section later.
          USE_FUNGICIDE   = .FALSE.
          FUNG_EFFICIENCY = 0.723
          FUNG_RES_D      = 14
          FUNG_BUF_D      = 16
          DVIP_THR        = 6
 
          CALL READ_DISEASE_PARAMETERS(CONTROL,
     &         LESION_S, KVERHULST, RVERHULST,
     &         YMAX, COF_A, COF_B,
     &         TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D, TMAX_D,
     &         LDMIN, LESIONAGEOPT, LESLIFEMAX, BETA, RRDS,
     &         NCYCLE, DAE_MIN, DAE_MIN_PRESENT,
     &         SRC_HALF, FAV_THR, NDS, DEP_FRAC)

          IS_MONOCYCLIC = (NCYCLE .EQ. 'M')

          SOURCE_PRESSURE    = 0.0
          SEC_SPORE_CLOUD    = 0.0
          SEC_SPORES_PENDING = 0.0
          FAV_SUM            = 0.0
          PRI_POOL           = 0.0
          PRI_RELEASED       = .FALSE.
          PEND_AREA          = 0.0
          PRESEASON_DONE     = .FALSE.
          PRESEASON_COUNT    = 0
          PRESEASON_DATE     = 0
          PRESEASON_RH       = 0.0
          PRESEASON_LWD      = 0.0
          PRESEASON_FT       = 0.0
          PRESEASON_FAV      = 0.0
 
          HDR_DONE        = .FALSE.
          OVERFLOW_WARNED = .FALSE.
 
          T   = 0.0
          LWD = 0.0
          FT  = 0.0
          FT_D = 0.0
          FT_G = 0.0
          DAILY_IP = 0.0
          HEALTH_LAI = 0.0
          LAI_SUSC   = 0.0
          DS_TOTAL   = 0.0
          LS_TODAY   = 0.0
          NEW_LOSS_TODAY = 0.0
          F_SOURCE = 0.0
          IS       = 0.0
 
          CALL DISMO_SEASON_RESET(MAXDAYS,
     &         ESP_LAT_HIST, ADMITTED_AREA, SLOT_DAE,
     &         DAE, PLANT_LIVE, N_ACTIVE_COH,
     &         LAI_PEAK_SEASON, CUM_NECROTIC, PREV_IS, SEVERITY_PCT,
     &         DVIP_pts, idx, SUM7, BufferDays, ResidualDays,
     &         NSprays, FungActive,
     &         DISEASE_LAI, DISEASE_SEN_RATE, VIRTUAL_PHOTO_FACTOR)
 
          CALL GETLUN('DISOUT', LUN_OUT)
          IF (CONTROL%RUN .EQ. 1) THEN
              OPEN(LUN_OUT, FILE='DISMO.OUT', STATUS='REPLACE',
     &             IOSTAT=IOS)
          ELSE
              OPEN(LUN_OUT, FILE='DISMO.OUT', STATUS='UNKNOWN',
     &             POSITION='APPEND', IOSTAT=IOS)
          ENDIF
          LUN_OPEN = (IOS .EQ. 0)
          IF (.NOT. LUN_OPEN) THEN
              MSG(1) = 'Cannot open DISMO.OUT. Disease output is'
              MSG(2) = 'disabled; the simulation continues.'
              CALL WARNING(2, 'DISMO ', MSG)
          ELSEIF (CONTROL%RUN .EQ. 1) THEN
              WRITE(LUN_OUT,'(A)')
     &         '*DISEASE IMPACT AND SEVERITY MODULE OUTPUT FILE'
          ENDIF
 
C***********************************************************************
C  SEASINIT -- once per season
C***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASINIT) THEN
 
          CALL DISMO_SEASON_RESET(MAXDAYS,
     &         ESP_LAT_HIST, ADMITTED_AREA, SLOT_DAE,
     &         DAE, PLANT_LIVE, N_ACTIVE_COH,
     &         LAI_PEAK_SEASON, CUM_NECROTIC, PREV_IS, SEVERITY_PCT,
     &         DVIP_pts, idx, SUM7, BufferDays, ResidualDays,
     &         NSprays, FungActive,
     &         DISEASE_LAI, DISEASE_SEN_RATE, VIRTUAL_PHOTO_FACTOR)
 
          SOURCE_PRESSURE    = 0.0
          SEC_SPORE_CLOUD    = 0.0
          SEC_SPORES_PENDING = 0.0
          FAV_SUM            = 0.0
          PRI_POOL           = 0.0
          PRI_RELEASED       = .FALSE.
          PEND_AREA          = 0.0
          PRESEASON_COUNT    = 0
          PRESEASON_DONE     = .FALSE.
          HDR_DONE           = .FALSE.
 
C  Re-open for seasons after the first in a sequential run.
          IF (.NOT. LUN_OPEN) THEN
              OPEN(LUN_OUT, FILE='DISMO.OUT', STATUS='UNKNOWN',
     &             POSITION='APPEND', IOSTAT=IOS)
              LUN_OPEN = (IOS .EQ. 0)
          END IF
 
          CALL DISMO_PRESEASON(CONTROL, TMIN_G, TOT_G, TMAX_G,
     &         TMIN_D, TOT_D, TMAX_D, USE_WTH_RH,
     &         SRC_SURV, SOURCE_PRESSURE, FAV_SUM,
     &         PRESEASON_DATE, PRESEASON_RH, PRESEASON_LWD,
     &         PRESEASON_FT, PRESEASON_FAV, PRESEASON_COUNT,
     &         MAXPRESEASON, INOC_LOOKBACK)
          PRESEASON_DONE = .TRUE.
 
C***********************************************************************
C  RATE 
C***********************************************************************
      ELSEIF (DYNAMIC .EQ. RATE) THEN
 
C-----------------------------------------------------------------------
C  1. ENVIRONMENT.  Unconditional, every day, canopy or no canopy.
C-----------------------------------------------------------------------

          CALL DISMO_UPDATE_ENVIRONMENT(Tmin, Tmax, RH, USE_WTH_RH,
     &         TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D, TMAX_D,
     &         SRC_SURV, DAILY_IP, SOURCE_PRESSURE, FAV_SUM,
     &         (CONTROL%DAS .GT. NVEG0),
     &         T, LWD, FT, FT_D, FT_G)
 
C  Airborne inoculum ages every day, canopy or no canopy.  Both
C  stores are aged in the same place so neither can freeze on a day
C  when deposition happens to be impossible.  Yesterday's emission
C  becomes airborne now.
          PRI_POOL        = PRI_POOL * SPOR_DECAY
          SEC_SPORE_CLOUD = SEC_SPORE_CLOUD * SEC_DECAY
     &                      + SEC_SPORES_PENDING
          SEC_SPORES_PENDING = 0.0
 
C-----------------------------------------------------------------------
C  2. EMERGENCE GATE AND DAE CLOCK
C-----------------------------------------------------------------------
          IF ((CONTROL%DAS .GT. NVEG0) .AND. (PLANT_LIVE .EQ. 0)) THEN
              DAE        = 1
              PLANT_LIVE = 1
          END IF
 
          IF (PLANT_LIVE .EQ. 1) THEN
 
C-----------------------------------------------------------------------
C  3. CANOPY BOOKKEEPING
C-----------------------------------------------------------------------
          LAI_PEAK_SEASON = MAX(LAI_PEAK_SEASON, LAI_TOTAL)

          LAI_SUSC   = MAX(LAI_PEAK_SEASON - CUM_NECROTIC, 0.0)
          IF (LAI_PEAK_SEASON .GT. EPS) THEN
              HEALTH_LAI = LAI_TOTAL
     &                   * MAX(1.0 - CUM_NECROTIC / LAI_PEAK_SEASON,0.0)
          ELSE
              HEALTH_LAI = LAI_TOTAL
          END IF
          HEALTH_LAI = MAX(HEALTH_LAI, 0.0)
 
          DAE_IDX       = MOD(DAE - 1, MAXDAYS) + 1
          EPI_ACTIVE    = (N_ACTIVE_COH .GT. 0)
          DEPOSITION_OK = (LAI_TOTAL .GE. LAI_MIN_START) .AND.
     &                    (HEALTH_LAI .GT. EPS)
 
C-----------------------------------------------------------------------
C  4. RISK INDEX AND FUNGICIDE.
C-----------------------------------------------------------------------
          CALL CALC_DVIP(LWD, T, DVIP_today)
          CALL APPLY_FUNGICIDE(DVIP_today, DVIP_pts, idx, SUM7,
     &         BufferDays, FungActive, ResidualDays, NSprays,
     &         USE_FUNGICIDE, FUNG_RES_D, FUNG_BUF_D, DVIP_THR,
     &         HEALTH_LAI)
 
C-----------------------------------------------------------------------
C  5. DEPOSITION AND NEW INFECTIONS
C-----------------------------------------------------------------------
          DS_TOTAL  = 0.0
          DS_PRI    = 0.0
          DS_SEC    = 0.0
          LS_TODAY  = 0.0
          FSS       = 0.0
          IR        = 0.0
          PRI_CLOUD = 0.0
          SEC_CLOUD = 0.0
          F_SOURCE  = 0.0

          IF (DEPOSITION_OK) THEN

              F_SOURCE = SOURCE_PRESSURE /
     &                   (SOURCE_PRESSURE + MAX(SRC_HALF, EPS))

C  --- primary inoculum arrival --------------------------------------

              IF (.NOT. PRI_RELEASED) THEN
                  IF (DAE_MIN_PRESENT) THEN
                      PRI_RELEASED = (DAE .GE. DAE_MIN)
                  ELSE
                      PRI_RELEASED = (FAV_SUM .GE. FAV_THR)
                  END IF
              END IF
              IF (PRI_RELEASED) THEN
                  PRI_POOL = PRI_POOL
     &                     + MAX(NDS * F_SOURCE * DAILY_IP, 0.0)
              END IF
              PRI_CLOUD = MAX(PRI_POOL, 0.0)

              SEC_CLOUD = MAX(SEC_SPORE_CLOUD, 0.0)
              CLOUD_TOTAL = PRI_CLOUD + SEC_CLOUD

              IF (CLOUD_TOTAL .GT. EPS) THEN
C  Both sources compete for the green area present today;
C  DS = MIN(cloud, canopy interception capacity).
                  CALL F_CANSPO(HEALTH_LAI, CLOUD_TOTAL, LESION_S,
     &                          DEP_FRAC, FSS)
                  CALL F_DS(FSS, CLOUD_TOTAL, DS_TOTAL)
                  DS_PRI = DS_TOTAL * PRI_CLOUD / CLOUD_TOTAL
                  DS_SEC = DS_TOTAL * SEC_CLOUD / CLOUD_TOTAL
C  Both pools are finite local quantities, so both are drawn down by
C  what lands on the canopy and by their own daily decay.
                  SEC_SPORE_CLOUD = MAX(SEC_SPORE_CLOUD - DS_SEC, 0.0)
                  PRI_POOL        = MAX(PRI_POOL        - DS_PRI, 0.0)
              END IF

              CALL F_IR(FT, LWD, YMAX, COF_A, COF_B, IR)
              IF (FungActive) IR = IR * (1.0 - FUNG_EFFICIENCY)
              CALL F_LS(IR, DS_TOTAL, LS_TODAY)
 
          END IF
 
C-----------------------------------------------------------------------
C  6. COHORT PROGRESSION.  Runs whenever cohorts exist or new
C     infections arrive
C-----------------------------------------------------------------------
          IS               = 0.0
          ESP_INOC_SEC     = 0.0
          NEW_LOSS_TODAY   = 0.0
          POT_SPO_PER_AREA = 0.0
          PEND_TODAY       = 0.0
          CLAIM_TODAY      = 0.0
 
          IF (EPI_ACTIVE .OR. (LS_TODAY .GT. 0.0)) THEN
 
C  Latency and ageing rates are constant within the day and do not
C  depend on k -- hoisted out of the loop 
              CALL F_LAR(FT_D, LESLIFEMAX, Lesion_Rate)
              CALL F_LR (FT_D, LDMIN, LR)
 
              IF (.NOT. IS_MONOCYCLIC) THEN
                  LES_COUNT_PREV = MAX(PREV_IS, 0.0)
     &                             / MAX(LESION_S, EPS)
                  CALL F_PPSR_POP(KVERHULST, RVERHULST,
     &                 LES_COUNT_PREV, PREV_IS, LAI_SUSC,
     &                 POT_SPO_PER_AREA)
              END IF
 
              DO k = 1, MAXDAYS
 
                  IF (SLOT_DAE(k) .LE. 0) CYCLE
 
                  ESP_LAT_HIST(k, C_LATP) = ESP_LAT_HIST(k, C_LATP)+LR
 
C  --- latent to infectious (once) ---
C  ADMITTED_AREA(k) is the cohort's FINAL footprint: the area those
C  lesions occupy once fully expanded.  It is fixed here, capped by
C  the leaf area still available, and drives sporulation from this day
C  on exactly as before.  It is NOT charged to the necrosis here.  
 
                  IF ((ESP_LAT_HIST(k, C_LATP) .GE. 1.0) .AND.
     &                (ESP_LAT_HIST(k, C_INFF) .LT. 0.5)) THEN

                      ESP_LAT_HIST(k, C_INFF) = 1.0
                      INF_AREA_K = MAX(ESP_LAT_HIST(k, C_LES), 0.0)
     &                             * LESION_S
                      LAI_AVAIL  = MAX(LAI_SUSC - PEND_AREA
     &                                 - CLAIM_TODAY, 0.0)
                      ADMITTED_AREA(k) = MIN(INF_AREA_K, LAI_AVAIL)
                      CLAIM_TODAY = CLAIM_TODAY + ADMITTED_AREA(k)
                  END IF
 
C  --- infectious phase ---
                  IF (ESP_LAT_HIST(k, C_INFF) .GT. 0.5) THEN
                      ESP_LAT_HIST(k, C_AGE) = ESP_LAT_HIST(k, C_AGE)
     &                                         + Lesion_Rate
                      LA = ESP_LAT_HIST(k, C_AGE)
 
                      IF (LA .GT. 1.0) THEN
C  Lifespan exhausted: retire the slot.  The necrotic area it created
C  stays in CUM_NECROTIC -- necrosis is permanent.
                          CALL DISMO_CLEAR_SLOT(MAXDAYS, k,
     &                         ESP_LAT_HIST, ADMITTED_AREA, SLOT_DAE)
                          N_ACTIVE_COH = MAX(N_ACTIVE_COH - 1, 0)
                          CYCLE
                      END IF
 
                      CALL F_LAF(LA, LAF, LESIONAGEOPT)
                      INF_AREA_K = ADMITTED_AREA(k)
                      IS = IS + INF_AREA_K

C  --- lesion expansion: this cohort's share of today's necrosis ---
C  Only the INCREMENT of the expanded fraction turns necrotic today.
C  Over the cohort's life the increments sum to ADMITTED_AREA(k)
C  exactly (expansion always completes before the lesion retires)
C  so the total necrosis a cohort produces is unchanged.  Only its
C  arrival is spread over the expansion period.
                      CALL F_LEXP(LA, LESIONAGEOPT, EXPN)
                      CALL F_LEXP(LA - Lesion_Rate, LESIONAGEOPT,
     &                            EXPN_PREV)
                      AREA_GROW = ADMITTED_AREA(k)
     &                            * MAX(EXPN - EXPN_PREV, 0.0)
                      AREA_GROW = MIN(AREA_GROW,
     &                            MAX(LAI_SUSC - NEW_LOSS_TODAY, 0.0))
                      NEW_LOSS_TODAY = NEW_LOSS_TODAY + AREA_GROW

C  What this cohort has claimed but not yet turned necrotic.  Summed
C  over live cohorts it is the area the canopy owes the epidemic, and
C  it is what the admission cap above must reserve.
                      PEND_TODAY = PEND_TODAY
     &                           + ADMITTED_AREA(k) * MAX(1.0-EXPN,0.0)
 
                      IF ((.NOT. IS_MONOCYCLIC) .AND.
     &                    (POT_SPO_PER_AREA .GT. 0.0) .AND.
     &                    (LAI_SUSC .GT. NEW_LOSS_TODAY)) THEN
                          PS_K = POT_SPO_PER_AREA * INF_AREA_K
     &                           * FT_D * LAF
                          ESP_INOC_SEC = ESP_INOC_SEC
     &                                   + MAX(PS_K, 0.0)
                      END IF
                  END IF
              END DO
 
C  --- register today's new infections ---
              IF (LS_TODAY .GT. 0.0) THEN
                  IF (SLOT_DAE(DAE_IDX) .GT. 0) THEN

                      IF (.NOT. OVERFLOW_WARNED) THEN
                          OVERFLOW_WARNED = .TRUE.
                          MSG(1) = 'Cohort ring buffer wrapped onto'
                          MSG(2) = 'a live cohort. Raise MAXDAYS or'
                          MSG(3) = 'lower LESLIFEMAX.'
                          CALL WARNING(3, 'DISMO ', MSG)
                      END IF
                      N_ACTIVE_COH = MAX(N_ACTIVE_COH - 1, 0)
                  END IF
                  CALL DISMO_CLEAR_SLOT(MAXDAYS, DAE_IDX,
     &                 ESP_LAT_HIST, ADMITTED_AREA, SLOT_DAE)
                  ESP_LAT_HIST(DAE_IDX, C_LES) = LS_TODAY
                  SLOT_DAE(DAE_IDX) = DAE
                  N_ACTIVE_COH = N_ACTIVE_COH + 1
              END IF
 
              IF (.NOT. IS_MONOCYCLIC) THEN
                  SEC_SPORES_PENDING = MAX(ESP_INOC_SEC, 0.0)
              ELSE
                  SEC_SPORES_PENDING = 0.0
              END IF
 
          END IF
 
C-----------------------------------------------------------------------
C  7. IMPACT.  necrosis, severity, photosynthetic
C     reduction and disease senescence are reported every day the crop
C     is alive.
C-----------------------------------------------------------------------
          CUM_NECROTIC = MIN(CUM_NECROTIC + NEW_LOSS_TODAY,
     &                       LAI_PEAK_SEASON)
          PREV_IS      = IS
          PEND_AREA    = PEND_TODAY
 
 
          CALL F_SEVERITY(CUM_NECROTIC, LAI_PEAK_SEASON, SEVERITY_PCT)
 
C  DISEASE_LAI exported in cm2 m-2
          IF (LAI_PEAK_SEASON .GT. EPS) THEN
              DISEASE_LAI = (CUM_NECROTIC / LAI_PEAK_SEASON)
     &                      * LAI_TOTAL * 10000.0
              DISEASE_LAI = MIN(DISEASE_LAI, LAI_TOTAL * 10000.0)
          ELSE
              DISEASE_LAI = 0.0
          END IF
          DISEASE_LAI = MAX(DISEASE_LAI, 0.0)
 
C  Virtual lesions
          SEVFRAC = MIN(MAX(SEVERITY_PCT / 100.0, 0.0), 1.0)
          CALL F_VIRTUAL_LESIONS(SEVFRAC, BETA, FVL)
          VIRTUAL_PHOTO_FACTOR = FVL
 
C  Disease senescence is driven by the actual severity.  
          CALL F_DEFOLIATION(WTLF, SLDOT, SEVERITY_PCT, RRDS,
     &                       DISEASE_SEN_RATE)
 
          ELSE
C  No live crop: hold the exported values neutral instead of stale.
              HEALTH_LAI           = LAI_TOTAL
              LAI_SUSC             = LAI_TOTAL
              DISEASE_LAI          = 0.0
              DISEASE_SEN_RATE     = 0.0
              VIRTUAL_PHOTO_FACTOR = 1.0
          END IF
 
          IF (PLANT_LIVE .EQ. 1) DAE = DAE + 1
          IF (YREND .EQ. YRDOY)  PLANT_LIVE = 0
 
C***********************************************************************
C  OUTPUT
C***********************************************************************
      ELSEIF (DYNAMIC .EQ. OUTPUT) THEN
 
        IF (LUN_OPEN) THEN
 
          IYEAR = YRDOY / 1000
          IDOY  = YRDOY - IYEAR * 1000
          IF (IYEAR .LT. 100) THEN
              IF (IYEAR .GT. 50) THEN
                  IYEAR = IYEAR + 1900
              ELSE
                  IYEAR = IYEAR + 2000
              END IF
          END IF
 
          IDAP = MAX(DAE - 1, 0)
 
          RH_OUT  = RH
          LWD_OUT = LWD
          FT_OUT  = FT
          FAV_OUT = FAV_SUM
          DO I_PRE = 1, PRESEASON_COUNT
              IF (PRESEASON_DATE(I_PRE) .EQ. YRDOY) THEN
                  RH_OUT  = PRESEASON_RH(I_PRE)
                  LWD_OUT = PRESEASON_LWD(I_PRE)
                  FT_OUT  = PRESEASON_FT(I_PRE)
                  FAV_OUT = PRESEASON_FAV(I_PRE)
                  EXIT
              END IF
          END DO
 
          IF (.NOT. HDR_DONE) THEN
              HDR_DONE = .TRUE.
              WRITE(LUN_OUT,'(A)') ' '
              WRITE(LUN_OUT,'(A,I4,A,A,A,A,1X,I5)')
     &         '*RUN',CONTROL%RUN,'        : ',
     &         CONTROL%ENAME,
     &         '                      ',
     &         CONTROL%MODEL,
     &         CONTROL%TRTNUM
              WRITE(LUN_OUT,'(A,A)')
     &         ' MODEL          : ', CONTROL%MODEL
              WRITE(LUN_OUT,'(A,A)')
     &         ' EXPERIMENT     : ', CONTROL%FILEX
              WRITE(LUN_OUT,'(A)')
     &         ' DATA PATH      :'
              WRITE(LUN_OUT,'(A,I2,A,A,A,A)')
     &         ' TREATMENT',CONTROL%TRTNUM,
     &         '    : ',CONTROL%ENAME,
     &         '                      ',CONTROL%MODEL
              WRITE(LUN_OUT,'(A)') ' '
              WRITE(LUN_OUT,'(A)') '!'
              WRITE(LUN_OUT,'(A)') '!'
              WRITE(LUN_OUT, 25)
          ENDIF
 
   25 FORMAT('@YEAR  DOY   DAS   DAE',
     &       '      LAIH      LWDh      RHU%      FTMP      LAIT',
     &       '      SUM7 NSPRAYS     FACT     SEV%',
     &       '        SRCP        FAVS        PRIM        SSCL',
     &       '        DSTO        LSTO        NLTO',
     &       '      FSRC')

          WRITE(LUN_OUT, 30)
     &     IYEAR, IDOY, CONTROL%DAS, IDAP,
     &     HEALTH_LAI,
     &     LWD_OUT, RH_OUT, FT_OUT, LAI_TOTAL, REAL(SUM7),
     &     NSprays, FungActive, SEVERITY_PCT,
     &     SOURCE_PRESSURE, FAV_OUT, PRI_CLOUD, SEC_SPORE_CLOUD,
     &     DS_TOTAL, LS_TODAY, NEW_LOSS_TODAY,
     &     F_SOURCE

   30 FORMAT(I5, I5, I6, I6, 5F10.3, F10.1, I8, L9, F9.2, 7F12.3,
     &       F10.4)
 
        ENDIF
 
C***********************************************************************
C  SEASEND
C***********************************************************************
      ELSEIF (DYNAMIC .EQ. SEASEND) THEN
 
          CALL DISMO_SEASON_RESET(MAXDAYS,
     &         ESP_LAT_HIST, ADMITTED_AREA, SLOT_DAE,
     &         DAE, PLANT_LIVE, N_ACTIVE_COH,
     &         LAI_PEAK_SEASON, CUM_NECROTIC, PREV_IS, SEVERITY_PCT,
     &         DVIP_pts, idx, SUM7, BufferDays, ResidualDays,
     &         NSprays, FungActive,
     &         DISEASE_LAI, DISEASE_SEN_RATE, VIRTUAL_PHOTO_FACTOR)
 
          IF (LUN_OPEN) THEN
              CLOSE(LUN_OUT)
              LUN_OPEN = .FALSE.
          END IF
 
      ENDIF
 
      RETURN
      END SUBROUTINE DISEASE_LEAF
 
!=======================================================================
!  SUBROUTINES
!=======================================================================
 
!-----------------------------------------------------------------------
!  SEASON RESET
!-----------------------------------------------------------------------
      SUBROUTINE DISMO_SEASON_RESET(MAXDAYS,
     &    ESP_LAT_HIST, ADMITTED_AREA, SLOT_DAE,
     &    DAE, PLANT_LIVE, N_ACTIVE_COH,
     &    LAI_PEAK_SEASON, CUM_NECROTIC, PREV_IS, SEVERITY_PCT,
     &    DVIP_pts, idx, SUM7, BufferDays, ResidualDays,
     &    NSprays, FungActive,
     &    DISEASE_LAI, DISEASE_SEN_RATE, VIRTUAL_PHOTO_FACTOR)
 
      IMPLICIT NONE
      INTEGER MAXDAYS
      REAL    ESP_LAT_HIST(MAXDAYS,4)
      REAL    ADMITTED_AREA(MAXDAYS)
      INTEGER SLOT_DAE(MAXDAYS)
      INTEGER DAE, PLANT_LIVE, N_ACTIVE_COH
      REAL    LAI_PEAK_SEASON, CUM_NECROTIC, PREV_IS, SEVERITY_PCT
      INTEGER DVIP_pts(7), idx, SUM7
      INTEGER BufferDays, ResidualDays, NSprays
      LOGICAL FungActive
      REAL    DISEASE_LAI, DISEASE_SEN_RATE, VIRTUAL_PHOTO_FACTOR
 
      ESP_LAT_HIST  = 0.0
      ADMITTED_AREA = 0.0
      SLOT_DAE      = 0
 
      DAE          = 0
      PLANT_LIVE   = 0
      N_ACTIVE_COH = 0
 
      LAI_PEAK_SEASON = 0.0
      CUM_NECROTIC    = 0.0
      PREV_IS         = 0.0
      SEVERITY_PCT    = 0.0
 
      DVIP_pts     = 0
      idx          = 1
      SUM7         = 0
      BufferDays   = 0
      ResidualDays = 0
      NSprays      = 0
      FungActive   = .FALSE.
 
      DISEASE_LAI          = 0.0
      DISEASE_SEN_RATE     = 0.0
      VIRTUAL_PHOTO_FACTOR = 1.0
 
      RETURN
      END SUBROUTINE DISMO_SEASON_RESET
 
!-----------------------------------------------------------------------
!  Clear one cohort slot of the ring buffer.
!-----------------------------------------------------------------------
      SUBROUTINE DISMO_CLEAR_SLOT(MAXDAYS, K,
     &                            ESP_LAT_HIST, ADMITTED_AREA, SLOT_DAE)
      IMPLICIT NONE
      INTEGER MAXDAYS, K
      REAL    ESP_LAT_HIST(MAXDAYS,4)
      REAL    ADMITTED_AREA(MAXDAYS)
      INTEGER SLOT_DAE(MAXDAYS)
 
      IF (K .LT. 1 .OR. K .GT. MAXDAYS) RETURN
      ESP_LAT_HIST(K,1) = 0.0
      ESP_LAT_HIST(K,2) = 0.0
      ESP_LAT_HIST(K,3) = 0.0
      ESP_LAT_HIST(K,4) = 0.0
      ADMITTED_AREA(K)  = 0.0
      SLOT_DAE(K)       = 0
 
      RETURN
      END SUBROUTINE DISMO_CLEAR_SLOT
!-----------------------------------------------------------------------
      SUBROUTINE DISMO_SPLIT(LINE, MAXTOK, TOKEN, NTOK)
      IMPLICIT NONE
      CHARACTER*(*)     LINE
      INTEGER           MAXTOK, NTOK
      CHARACTER(LEN=32) TOKEN(MAXTOK)
      INTEGER I, J, L
      LOGICAL INTOK
      CHARACTER(LEN=1) C
 
      NTOK  = 0
      J     = 0
      INTOK = .FALSE.
      DO I = 1, MAXTOK
          TOKEN(I) = ' '
      END DO
 
      L = LEN_TRIM(LINE)
      DO I = 1, L
          C = LINE(I:I)
          IF (C .EQ. ' ' .OR. C .EQ. CHAR(9) .OR. C .EQ. CHAR(13)) THEN
              INTOK = .FALSE.
          ELSE
              IF (.NOT. INTOK) THEN
                  IF (NTOK .GE. MAXTOK) RETURN
                  NTOK  = NTOK + 1
                  J     = 0
                  INTOK = .TRUE.
              END IF
              IF (J .LT. 32) THEN
                  J = J + 1
                  TOKEN(NTOK)(J:J) = C
              END IF
          END IF
      END DO
 
      RETURN
      END SUBROUTINE DISMO_SPLIT
 
!-----------------------------------------------------------------------
!  Read a REAL from token IPOS.  Missing token or -99 -> DEFVAL.
!-----------------------------------------------------------------------
      SUBROUTINE DISMO_TOKR(TOKEN, MAXTOK, NTOK, IPOS, DEFVAL, VALUE,
     &                      OK)
      IMPLICIT NONE
      INTEGER           MAXTOK, NTOK, IPOS
      CHARACTER(LEN=32) TOKEN(MAXTOK)
      REAL              DEFVAL, VALUE
      LOGICAL           OK
      INTEGER           IOS
      REAL              TMP
 
      VALUE = DEFVAL
      OK    = .FALSE.
      IF (IPOS .LT. 1 .OR. IPOS .GT. NTOK) RETURN
      IF (LEN_TRIM(TOKEN(IPOS)) .EQ. 0) RETURN
      READ(TOKEN(IPOS), *, IOSTAT=IOS) TMP
      IF (IOS .NE. 0) RETURN
      IF (ABS(TMP + 99.0) .LT. 1.0E-4) RETURN
      VALUE = TMP
      OK    = .TRUE.
 
      RETURN
      END SUBROUTINE DISMO_TOKR
 
!-----------------------------------------------------------------------
!  READ DISEASE PARAMETERS
!
!  DAE_START = -99 removes the window guard: onset is then decided
!  purely by accumulated weather.  A blank field is NOT accepted --
!  the reader is whitespace based, so a blank would shift every
!  following column by one position.
!-----------------------------------------------------------------------
      SUBROUTINE READ_DISEASE_PARAMETERS(CONTROL,
     &    LESION_S, KVERHULST, RVERHULST,
     &    YMAX, COF_A, COF_B,
     &    TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D, TMAX_D,
     &    LDMIN, LESIONAGEOPT, LESLIFEMAX, BETA, RRDS,
     &    NCYCLE, DAE_MIN, DAE_MIN_PRESENT,
     &    SRC_HALF, FAV_THR, NDS, DEP_FRAC)

      USE ModuleDefs
      IMPLICIT NONE
      EXTERNAL GETLUN, ERROR, WARNING, DISMO_SPLIT, DISMO_TOKR

      TYPE (ControlType) CONTROL

      REAL    LESION_S, KVERHULST, RVERHULST
      REAL    YMAX, COF_A, COF_B
      REAL    TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D, TMAX_D
      REAL    LDMIN, LESIONAGEOPT, LESLIFEMAX, BETA, RRDS
      REAL    SRC_HALF, FAV_THR, NDS, DEP_FRAC
      INTEGER DAE_MIN
      LOGICAL DAE_MIN_PRESENT
      CHARACTER(LEN=1) NCYCLE

      INTEGER, PARAMETER :: MAXTOK = 30
      INTEGER, PARAMETER :: NREQ   = 25
      CHARACTER(LEN=6),  PARAMETER :: ERRKEY = 'DISMO '
 
      CHARACTER(LEN=32)  TOKEN(MAXTOK)
      CHARACTER(LEN=400) LINE
      CHARACTER(LEN=32)  TARGET_DISEASE
      CHARACTER(LEN=120) DISFIL
      CHARACTER(LEN=78)  MSG(4)
      INTEGER LUN_DIS, IOS, STATE, NTOK, LNUM
      LOGICAL FEXIST, FOUND_DISEASE, OK
      REAL    RTMP
 
      DISFIL = 'disease_parameters.txt'
      INQUIRE(FILE=TRIM(DISFIL), EXIST=FEXIST)
      IF (.NOT. FEXIST) THEN
          MSG(1) = 'disease_parameters.txt was not found in the'
          MSG(2) = 'working directory.'
          CALL WARNING(2, ERRKEY, MSG)
          CALL ERROR(ERRKEY, 29, DISFIL, 0)
      END IF
 
      CALL GETLUN('DISINP', LUN_DIS)
      OPEN(LUN_DIS, FILE=TRIM(DISFIL), STATUS='OLD', ACTION='READ',
     &     IOSTAT=IOS)
      IF (IOS .NE. 0) CALL ERROR(ERRKEY, IOS, DISFIL, 0)
 
!     STATE 0 : looking for a section tag
!     STATE 1 : inside *DISEASE CONTROL, waiting for the @ header
!     STATE 2 : next valid record holds the target disease name
!     STATE 3 : inside *DISEASE DATABASE, waiting for the @ header
!     STATE 4 : reading database records
      TARGET_DISEASE = ' '
      STATE          = 0
      FOUND_DISEASE  = .FALSE.
      DAE_MIN_PRESENT = .FALSE.
      DAE_MIN        = -99
      NCYCLE         = 'P'
      LNUM           = 0
      FAV_THR        = -99.0
      NDS             = 0.0
      DEP_FRAC        = 0.25
 
      DO WHILE (.TRUE.)
          READ(LUN_DIS, '(A)', IOSTAT=IOS) LINE
          IF (IOS .NE. 0) EXIT
          LNUM = LNUM + 1
          LINE = ADJUSTL(LINE)
          IF (LEN_TRIM(LINE) .EQ. 0) CYCLE
          IF (LINE(1:1) .EQ. '!') CYCLE
 
          SELECT CASE (STATE)
 
          CASE (0)
              IF (INDEX(LINE,'*DISEASE CONTROL') .GT. 0) THEN
                  STATE = 1
              ELSEIF (INDEX(LINE,'*DISEASE DATABASE') .GT. 0) THEN
                  STATE = 3
              END IF
 
          CASE (1)
              IF (LINE(1:1) .EQ. '@') STATE = 2
 
          CASE (2)
              CALL DISMO_SPLIT(LINE, MAXTOK, TOKEN, NTOK)
              IF (NTOK .GE. 1) TARGET_DISEASE = TOKEN(1)
              STATE = 0
 
          CASE (3)
              IF (LINE(1:1) .EQ. '@') STATE = 4
 
          CASE (4)
              CALL DISMO_SPLIT(LINE, MAXTOK, TOKEN, NTOK)
              IF (NTOK .LT. 2) CYCLE
              IF (TRIM(TOKEN(2)) .NE. TRIM(TARGET_DISEASE)) CYCLE
 
!             A record that fills the buffer was probably truncated,
!             which would silently drop the trailing optional columns
!             and replace them with defaults.
              IF (LEN_TRIM(LINE) .GE. 400) THEN
                  MSG(1) = 'Disease record reaches the 400-character'
                  MSG(2) = 'buffer and may be truncated. Trailing'
                  MSG(3) = 'trailing columns may be lost.'
                  CALL WARNING(3, ERRKEY, MSG)
              END IF
 
              IF (NTOK .LT. NREQ) THEN
                  MSG(1) = 'Disease record has fewer than 25 required'
                  MSG(2) = 'columns. CLD_THR and CLD_W were removed,'
                  MSG(3) = 'FAV_THR added and NDS redefined: rewrite'
                  MSG(4) = 'the file, do not patch it.'
                  CALL WARNING(4, ERRKEY, MSG)
                  CALL ERROR(ERRKEY, 59, DISFIL, LNUM)
              END IF

!             Column 3 is the observed onset date in days after
!             emergence, or -99 when there is none and onset is left to
!             the weather.  DISMO_TOKR reports -99 as "not present".
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK, 3, -99.0, RTMP,  OK)
              IF (OK) THEN
                  DAE_MIN         = NINT(RTMP)
                  DAE_MIN_PRESENT = .TRUE.
              END IF
              NCYCLE = TOKEN(4)(1:1)
              IF (NCYCLE .NE. 'M' .AND. NCYCLE .NE. 'P') NCYCLE = 'P'

              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK, 5,  1.0, SRC_HALF,OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK, 6, 30.0, FAV_THR, OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK, 7,  0.0, NDS,   OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK, 8, 0.25, DEP_FRAC,
     &                        OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK, 9,1.0E-6,LESION_S,OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,10,  1.0, KVERHULST,
     &                        OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,11,  0.1, RVERHULST,
     &                        OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,12,  1.0, YMAX,  OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,13,  0.1, COF_A, OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,14,  1.0, COF_B, OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,15, 10.0, TMIN_G,OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,16, 22.0, TOT_G, OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,17, 30.0, TMAX_G,OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,18, 10.0, TMIN_D,OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,19, 24.0, TOT_D, OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,20, 32.0, TMAX_D,OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,21,  6.0, LDMIN, OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,22,  0.5, LESIONAGEOPT,
     &                        OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,23, 30.0, LESLIFEMAX,
     &                        OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,24,  1.0, BETA,  OK)
              CALL DISMO_TOKR(TOKEN,MAXTOK,NTOK,25,  0.0, RRDS,  OK)

              FOUND_DISEASE = .TRUE.
              EXIT
 
          END SELECT
      END DO
 
      CLOSE(LUN_DIS)
 
      IF (STATE .LT. 4) THEN
          MSG(1) = 'No *DISEASE DATABASE section with an @ header was'
          MSG(2) = 'found in disease_parameters.txt.'
          CALL WARNING(2, ERRKEY, MSG)
          CALL ERROR(ERRKEY, 59, DISFIL, LNUM)
      END IF
 
      IF (.NOT. FOUND_DISEASE) THEN
          MSG(1) = 'Target disease not found in the database:'
          MSG(2) = TRIM(TARGET_DISEASE)
          MSG(3) = 'Check the *DISEASE CONTROL section.'
          CALL WARNING(3, ERRKEY, MSG)
          CALL ERROR(ERRKEY, 59, DISFIL, LNUM)
      END IF
 
!     Basic sanity checks: these catch parameter sets that would make
!     the temperature or lesion-age responses undefined.
      IF (TOT_G .LE. TMIN_G .OR. TMAX_G .LE. TOT_G .OR.
     &    TOT_D .LE. TMIN_D .OR. TMAX_D .LE. TOT_D) THEN
          MSG(1) = 'Cardinal temperatures must satisfy'
          MSG(2) = 'TMIN < TOPT < TMAX for both G and D responses.'
          CALL WARNING(2, ERRKEY, MSG)
          CALL ERROR(ERRKEY, 59, DISFIL, LNUM)
      END IF
 
      IF (LESIONAGEOPT .LE. 0.0 .OR. LESIONAGEOPT .GE. 1.0) THEN
          MSG(1) = 'LESIONAGEOPT must be a fraction of LESLIFEMAX in'
          MSG(2) = 'the open interval (0,1). Value clamped to 0.5.'
          CALL WARNING(2, ERRKEY, MSG)
          LESIONAGEOPT = 0.5
      END IF
 
      IF (BETA .LT. 1.0) THEN
          MSG(1) = 'BETA below 1 implies a virtual lesion smaller than'
          MSG(2) = 'the visible lesion. Clamped to 1 (no extra effect).'
          CALL WARNING(2, ERRKEY, MSG)
          BETA = 1.0
      END IF

!     Without a dose the arrival event releases nothing and the epidemic
!     never starts, which reads as a model failure rather than as the
!     misconfigured file it is.
      IF (NDS .LE. 0.0) THEN
          MSG(1) = 'NDS must be positive: with no arriving spores the'
          MSG(2) = 'inoculum arrival releases nothing and no epidemic'
          MSG(3) = 'can start. Check column 7.'
          CALL WARNING(3, ERRKEY, MSG)
          CALL ERROR(ERRKEY, 59, DISFIL, LNUM)
      END IF

!     A settling fraction outside (0,1] either freezes the cloud or
!     deposits more than exists.
      IF (DEP_FRAC .LE. 0.0 .OR. DEP_FRAC .GT. 1.0) THEN
          MSG(1) = 'DEP_FRAC must lie in (0,1]: it is the fraction'
          MSG(2) = 'of the airborne cloud that settles in one day.'
          MSG(3) = 'Use 1.0 for the old instantaneous deposition.'
          CALL WARNING(3, ERRKEY, MSG)
          CALL ERROR(ERRKEY, 59, DISFIL, LNUM)
      END IF

!     FAV_THR is only read when there is no observed onset date.
      IF (.NOT. DAE_MIN_PRESENT .AND. FAV_THR .LE. 0.0) THEN
          MSG(1) = 'DAE_START is -99, so onset is weather driven and'
          MSG(2) = 'FAV_THR must be positive. Set an observed onset in'
          MSG(3) = 'DAE_START or a favourability threshold in FAV_THR.'
          CALL WARNING(3, ERRKEY, MSG)
          CALL ERROR(ERRKEY, 59, DISFIL, LNUM)
      END IF

      RETURN
      END SUBROUTINE READ_DISEASE_PARAMETERS
 
!-----------------------------------------------------------------------
!  PRE-PLANT ENVIRONMENTAL INOCULUM RECONSTRUCTION
!-----------------------------------------------------------------------
      SUBROUTINE DISMO_PRESEASON(CONTROL, TMIN_G, TOT_G, TMAX_G,
     &    TMIN_D, TOT_D, TMAX_D, USE_WTH_RH, SRC_SURV,
     &    SOURCE_PRESSURE, FAV_SUM,
     &    PRESEASON_DATE, PRESEASON_RH, PRESEASON_LWD,
     &    PRESEASON_FT, PRESEASON_FAV, PRESEASON_COUNT,
     &    MAXPRESEASON, INOC_LOOKBACK)

      USE ModuleDefs
      IMPLICIT NONE
      EXTERNAL GETLUN, WARNING, DISMO_UPDATE_ENVIRONMENT,
     &         DISMO_NORMALIZE_DATE, DISMO_READ_WTH_FIELD,
     &         DISMO_SHIFT_DATE

      TYPE (ControlType) CONTROL

      REAL    TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D, TMAX_D
      REAL    SRC_SURV, SOURCE_PRESSURE, FAV_SUM
      INTEGER MAXPRESEASON, PRESEASON_COUNT, INOC_LOOKBACK
      INTEGER PRESEASON_DATE(MAXPRESEASON)
      REAL    PRESEASON_RH(MAXPRESEASON)
      REAL    PRESEASON_LWD(MAXPRESEASON)
      REAL    PRESEASON_FT(MAXPRESEASON)
      REAL    PRESEASON_FAV(MAXPRESEASON)
      LOGICAL USE_WTH_RH
 
      CHARACTER(LEN=6), PARAMETER :: ERRKEY = 'DISMO '
 
      REAL    WTMAX, WTMIN, WRH, DAILY_IP, T, LWD, FT, FT_D, FT_G
      LOGICAL FOUND_HEADER, IN_PLANTING, OK_TMAX, OK_TMIN, OK_RH
      LOGICAL SEP_SLASH
      INTEGER LUN_IO, LUN_WTH, IOS, WTH_DATE, FULL_DATE
      INTEGER POS_TMAX, POS_TMIN, POS_RAIN, POS_RHUM
      INTEGER PDATE, START_DATE, NDAYS_READ
      CHARACTER(LEN=400) LINE, WTH_HEADER
      CHARACTER(LEN=30)  WTH_NAME
      CHARACTER(LEN=120) WTH_PATH
      CHARACTER(LEN=240) WTH_FILE
      CHARACTER(LEN=78)  MSG(4)
 
      PDATE       = -99
      WTH_NAME    = ' '
      WTH_PATH    = ' '
      IN_PLANTING = .FALSE.
      NDAYS_READ  = 0
 
!---- Resolved planting date and weather location ---------------------
      CALL GETLUN('DSMINP', LUN_IO)
      OPEN(LUN_IO, FILE=TRIM(CONTROL%FILEIO), STATUS='OLD',
     &     ACTION='READ', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
          MSG(1) = 'Cannot open the DSSAT input file.'
          MSG(2) = 'Pre-plant inoculum reconstruction skipped.'
          CALL WARNING(2, ERRKEY, MSG)
          RETURN
      END IF
 
      DO WHILE (.TRUE.)
          READ(LUN_IO, '(A)', IOSTAT=IOS) LINE
          IF (IOS .NE. 0) EXIT
          LINE = ADJUSTL(LINE)
          IF (LEN_TRIM(LINE) .EQ. 0) CYCLE
 
          IF (INDEX(LINE, 'WEATHERW') .EQ. 1) THEN
              READ(LINE(9:), *, IOSTAT=IOS) WTH_NAME, WTH_PATH
              CYCLE
          END IF
 
          IF (INDEX(LINE, '*PLANTING DETAILS') .EQ. 1) THEN
              IN_PLANTING = .TRUE.
              CYCLE
          END IF
 
          IF (IN_PLANTING) THEN
              IF (LINE(1:1) .EQ. '*') THEN
                  IN_PLANTING = .FALSE.
              ELSE
                  READ(LINE, *, IOSTAT=IOS) PDATE
                  IF (IOS .EQ. 0 .AND. PDATE .GT. 0) IN_PLANTING =
     &                .FALSE.
              END IF
          END IF
      END DO
      CLOSE(LUN_IO)
 
      IF (LEN_TRIM(WTH_NAME) .EQ. 0 .OR. PDATE .LE. 0) THEN
          MSG(1) = 'Weather file or planting date not resolved.'
          MSG(2) = 'Pre-plant inoculum reconstruction skipped.'
          CALL WARNING(2, ERRKEY, MSG)
          RETURN
      END IF
 
!     The replay window is anchored on planting, not on YRSIM.
      CALL DISMO_SHIFT_DATE(PDATE, INOC_LOOKBACK, START_DATE)
 
!---- Build the weather file path --------------------------------------
      SEP_SLASH = (INDEX(WTH_PATH, '/') .GT. 0)
      IF (LEN_TRIM(WTH_PATH) .EQ. 0) THEN
          WTH_FILE = TRIM(WTH_NAME)
      ELSEIF (WTH_PATH(LEN_TRIM(WTH_PATH):LEN_TRIM(WTH_PATH))
     &        .EQ. CHAR(92) .OR.
     &        WTH_PATH(LEN_TRIM(WTH_PATH):LEN_TRIM(WTH_PATH))
     &        .EQ. '/') THEN
          WTH_FILE = TRIM(WTH_PATH)//TRIM(WTH_NAME)
      ELSEIF (SEP_SLASH) THEN
          WTH_FILE = TRIM(WTH_PATH)//'/'//TRIM(WTH_NAME)
      ELSE
          WTH_FILE = TRIM(WTH_PATH)//CHAR(92)//TRIM(WTH_NAME)
      END IF
 
      CALL GETLUN('DSMWTH', LUN_WTH)
      OPEN(LUN_WTH, FILE=TRIM(WTH_FILE), STATUS='OLD', ACTION='READ',
     &     IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
          MSG(1) = 'Cannot open the weather file:'
          MSG(2) = TRIM(WTH_FILE)
          MSG(3) = 'Pre-plant inoculum reconstruction skipped.'
          CALL WARNING(3, ERRKEY, MSG)
          RETURN
      END IF
 
      FOUND_HEADER = .FALSE.
      POS_TMAX = 0
      POS_TMIN = 0
      POS_RAIN = 0
      POS_RHUM = 0
 
      DO WHILE (.TRUE.)
          READ(LUN_WTH, '(A)', IOSTAT=IOS) LINE
          IF (IOS .NE. 0) EXIT
          LINE = ADJUSTL(LINE)
 
          IF (.NOT. FOUND_HEADER) THEN
              IF (INDEX(LINE, '@DATE') .EQ. 1) THEN
                  WTH_HEADER = LINE
                  POS_TMAX = INDEX(WTH_HEADER, 'TMAX')
                  POS_TMIN = INDEX(WTH_HEADER, 'TMIN')
                  POS_RAIN = INDEX(WTH_HEADER, 'RAIN')
                  POS_RHUM = INDEX(WTH_HEADER, 'RHUM')
                  IF (POS_TMAX .GT. 0 .AND. POS_TMIN .GT. 0
     &                .AND. POS_RAIN .GT. 0) THEN
                      IF (USE_WTH_RH .AND. POS_RHUM .EQ. 0) THEN
                          MSG(1) = 'RHUM is missing from the weather'
                          MSG(2) = 'file. Pre-plant reconstruction'
                          MSG(3) = 'skipped.'
                          CALL WARNING(3, ERRKEY, MSG)
                          CLOSE(LUN_WTH)
                          RETURN
                      END IF
                      FOUND_HEADER = .TRUE.
                  END IF
              END IF
              CYCLE
          END IF
 
          IF (LEN_TRIM(LINE) .EQ. 0 .OR. LINE(1:1) .EQ. '!') CYCLE
          READ(LINE(1:5), '(I5)', IOSTAT=IOS) WTH_DATE
          IF (IOS .NE. 0) CYCLE
 
          CALL DISMO_NORMALIZE_DATE(WTH_DATE, FULL_DATE)
          IF (FULL_DATE .LT. START_DATE .OR. FULL_DATE .GE. PDATE)
     &        CYCLE
 
          CALL DISMO_READ_WTH_FIELD(LINE, POS_TMAX, POS_TMIN-1,
     &         WTMAX, OK_TMAX)
          CALL DISMO_READ_WTH_FIELD(LINE, POS_TMIN, POS_RAIN-1,
     &         WTMIN, OK_TMIN)
          IF (USE_WTH_RH) THEN
              CALL DISMO_READ_WTH_FIELD(LINE, POS_RHUM, LEN(LINE),
     &             WRH, OK_RH)
          ELSE
              WRH   = 0.0
              OK_RH = .TRUE.
          END IF
 
          IF (OK_TMAX .AND. OK_TMIN .AND. OK_RH) THEN
C  .FALSE.: the replay builds SOURCE_PRESSURE, which is its whole
C  purpose, but it must NOT advance the arrival clock -- there is no
C  crop to expose yet.  PRESEASON_FAV therefore stays zero, and the
C  FAVS column reads zero over the replayed days, which is honest.
              CALL DISMO_UPDATE_ENVIRONMENT(WTMIN, WTMAX, WRH,
     &             USE_WTH_RH, TMIN_G, TOT_G, TMAX_G,
     &             TMIN_D, TOT_D, TMAX_D, SRC_SURV,
     &             DAILY_IP, SOURCE_PRESSURE, FAV_SUM,
     &             .FALSE.,
     &             T, LWD, FT, FT_D, FT_G)
              NDAYS_READ = NDAYS_READ + 1

              IF (PRESEASON_COUNT .LT. MAXPRESEASON) THEN
                  PRESEASON_COUNT = PRESEASON_COUNT + 1
                  PRESEASON_DATE(PRESEASON_COUNT) = FULL_DATE
                  PRESEASON_RH(PRESEASON_COUNT)   = WRH
                  PRESEASON_LWD(PRESEASON_COUNT)  = LWD
                  PRESEASON_FT(PRESEASON_COUNT)   = FT
                  PRESEASON_FAV(PRESEASON_COUNT)  = FAV_SUM
              END IF
          END IF
      END DO
      CLOSE(LUN_WTH)
 
      IF (.NOT. FOUND_HEADER) THEN
          MSG(1) = 'Cannot read the @DATE header from:'
          MSG(2) = TRIM(WTH_FILE)
          CALL WARNING(2, ERRKEY, MSG)
      ELSEIF (NDAYS_READ .LT. INOC_LOOKBACK) THEN

          WRITE(MSG(1),'(A,I4,A,I4,A)')
     &     'Pre-plant replay covered ', NDAYS_READ, ' of ',
     &     INOC_LOOKBACK, ' days.'
          MSG(2) = 'Weather record starts too late; initial inoculum'
          MSG(3) = 'pressure is underestimated.'
          CALL WARNING(3, ERRKEY, MSG)
      END IF
 
      RETURN
      END SUBROUTINE DISMO_PRESEASON
 
!-----------------------------------------------------------------------
!  DAILY ENVIRONMENTAL UPDATE
!-----------------------------------------------------------------------
      SUBROUTINE DISMO_UPDATE_ENVIRONMENT(TMIN, TMAX, RH, USE_WTH_RH,
     &    TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D, TMAX_D,
     &    SRC_SURV, DAILY_IP, SOURCE_PRESSURE, FAV_SUM,
     &    COUNT_EXPOSURE, T, LWD, FT, FT_D, FT_G)

      IMPLICIT NONE
      EXTERNAL F_TAVG, F_DEW, F_RH, F_LWD, T_DEV

      REAL TMIN, TMAX, RH, TMIN_G, TOT_G, TMAX_G
      REAL TMIN_D, TOT_D, TMAX_D, SRC_SURV
      REAL DAILY_IP, SOURCE_PRESSURE, T, LWD, FT, FT_D
      REAL FAV_SUM
      REAL FT_G, TDEW, ES, E, RH_LOCAL, T_WET
      REAL FT_SRC, FT_D_SRC, FT_G_SRC
      LOGICAL USE_WTH_RH, COUNT_EXPOSURE

      CALL F_TAVG(TMAX, TMIN, T)
      RH_LOCAL = RH
      IF (.NOT. USE_WTH_RH) THEN
          CALL F_DEW(T, TMIN, TMAX, TDEW)
          CALL F_RH(RH_LOCAL, TDEW, ES, E, T)
      END IF
      CALL F_LWD(RH_LOCAL, LWD)

      T_WET = TMIN
      CALL T_DEV(T_WET, TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D, TMAX_D,
     &           FT, FT_D, FT_G)
 

      CALL T_DEV(T, TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D, TMAX_D,
     &           FT_SRC, FT_D_SRC, FT_G_SRC)
      DAILY_IP = MAX(FT_G_SRC, 0.0) * MIN(MAX(LWD / 24.0, 0.0), 1.0)
 
      SOURCE_PRESSURE = SOURCE_PRESSURE * SRC_SURV + DAILY_IP
      IF (COUNT_EXPOSURE) THEN
          FAV_SUM     = FAV_SUM                    + DAILY_IP
      END IF

      RETURN
      END SUBROUTINE DISMO_UPDATE_ENVIRONMENT
 
!-----------------------------------------------------------------------
!  Normalise YYDDD or YYYYDDD to YYYYDDD.
!-----------------------------------------------------------------------
      SUBROUTINE DISMO_NORMALIZE_DATE(INPUT_DATE, FULL_DATE)
      IMPLICIT NONE
      INTEGER INPUT_DATE, FULL_DATE, IYEAR, IDOY
 
      IYEAR = INPUT_DATE / 1000
      IDOY  = INPUT_DATE - IYEAR * 1000
      IF (IYEAR .LT. 100) THEN
          IF (IYEAR .GT. 50) THEN
              IYEAR = IYEAR + 1900
          ELSE
              IYEAR = IYEAR + 2000
          END IF
      END IF
      FULL_DATE = IYEAR * 1000 + IDOY
 
      RETURN
      END SUBROUTINE DISMO_NORMALIZE_DATE
 
!-----------------------------------------------------------------------
!  Shift a YYYYDDD date backward by OFFSET_DAYS.
!-----------------------------------------------------------------------
      SUBROUTINE DISMO_SHIFT_DATE(FULL_DATE, OFFSET_DAYS, NEW_DATE)
      IMPLICIT NONE
      INTEGER FULL_DATE, OFFSET_DAYS, NEW_DATE
      INTEGER IYEAR, IDOY, DAYS_IN_YEAR, DAYS_LEFT
      LOGICAL LEAP_YEAR
 
      IYEAR     = FULL_DATE / 1000
      IDOY      = FULL_DATE - IYEAR * 1000
      DAYS_LEFT = MAX(OFFSET_DAYS, 0)
 
      DO WHILE (DAYS_LEFT .GT. 0)
          IF (IDOY .GT. DAYS_LEFT) THEN
              IDOY      = IDOY - DAYS_LEFT
              DAYS_LEFT = 0
          ELSE
              DAYS_LEFT = DAYS_LEFT - IDOY
              IYEAR     = IYEAR - 1
              LEAP_YEAR = MOD(IYEAR,400) .EQ. 0 .OR.
     &          (MOD(IYEAR,4) .EQ. 0 .AND. MOD(IYEAR,100) .NE. 0)
              IF (LEAP_YEAR) THEN
                  DAYS_IN_YEAR = 366
              ELSE
                  DAYS_IN_YEAR = 365
              END IF
              IDOY = DAYS_IN_YEAR
          END IF
      END DO
 
      NEW_DATE = IYEAR * 1000 + IDOY
 
      RETURN
      END SUBROUTINE DISMO_SHIFT_DATE
 
!-----------------------------------------------------------------------
!  Read one fixed-width value from a DSSAT WTH record.
!-----------------------------------------------------------------------
      SUBROUTINE DISMO_READ_WTH_FIELD(LINE, FIRST, LAST, VALUE, OK)
      IMPLICIT NONE
      CHARACTER*(*) LINE
      INTEGER FIRST, LAST, LAST_POS, IOS
      REAL    VALUE
      LOGICAL OK
 
      OK    = .FALSE.
      VALUE = 0.0
      IF (FIRST .LT. 1 .OR. LAST .LT. FIRST) RETURN
      IF (FIRST .GT. LEN(LINE)) RETURN
      LAST_POS = MIN(LAST, LEN(LINE))
      READ(LINE(FIRST:LAST_POS), *, IOSTAT=IOS) VALUE
      IF (IOS .EQ. 0) OK = .TRUE.
 
      RETURN
      END SUBROUTINE DISMO_READ_WTH_FIELD
 
!=======================================================================
!  PROCESS FUNCTIONS
!=======================================================================
 
!----- Mean temperature ------------------------------------------------
      SUBROUTINE F_TAVG(Tmax, Tmin, T)
      IMPLICIT NONE
      REAL Tmax, Tmin, T
      T = (Tmax + Tmin) / 2.0
      RETURN
      END SUBROUTINE F_TAVG
 
!----- Dew point (empirical) -------------------------------------------
      SUBROUTINE F_DEW(T, Tmin, Tmax, Tdew)
      IMPLICIT NONE
      REAL T, Tmin, Tmax, Tdew
      Tdew = (-0.036*T) + (0.9679*Tmin) + (0.0072*(Tmax-Tmin)) + 1.0111
      RETURN
      END SUBROUTINE F_DEW
 
!----- Relative humidity from dew point (Tetens) -----------------------
      SUBROUTINE F_RH(RH, Tdew, Es, E, T)
      IMPLICIT NONE
      REAL RH, Tdew, Es, E, T
      Es = EXP(17.625*T    / (243.04 + T))
      E  = EXP(17.625*Tdew / (243.04 + Tdew))
      RH = MIN(MAX((E / Es) * 100.0, 0.0), 100.0)
      RETURN
      END SUBROUTINE F_RH
 
!----- Leaf wetness duration from RH -----------------------------------
      SUBROUTINE F_LWD(RH, LWD)
      IMPLICIT NONE
      REAL RH, LWD, RH_LOC
 
      RH_LOC = RH
      IF (RH_LOC .LE. 1.0) RH_LOC = RH_LOC * 100.0
      RH_LOC = MIN(MAX(RH_LOC, 0.0), 100.0)
      LWD = 31.31 / (1.0 + EXP(-((RH_LOC - 85.17) / 9.13)))
      LWD = MIN(MAX(LWD, 0.0), 24.0)
 
      RETURN
      END SUBROUTINE F_LWD
 
!----- Beta temperature responses --------------------------------------
      SUBROUTINE T_DEV(T, TMIN_G, TOT_G, TMAX_G, TMIN_D, TOT_D, TMAX_D,
     &                 FT, FT_D, FT_G)
      IMPLICIT NONE
      REAL FT_D, FT_G, FT, T
      REAL TMIN_D, TOT_D, TMAX_D, TMIN_G, TOT_G, TMAX_G
 
      FT_G = 0.0
      IF ((TMAX_G - T) .GT. 0.0 .AND. (T - TMIN_G) .GT. 0.0 .AND.
     &    (TMAX_G - TOT_G) .GT. 0.0 .AND. (TOT_G - TMIN_G) .GT. 0.0)
     &    THEN
          FT_G = ((TMAX_G-T)/(TMAX_G-TOT_G)) *
     &    ((T-TMIN_G)/(TOT_G-TMIN_G))**((TOT_G-TMIN_G)/(TMAX_G-TOT_G))
      END IF
 
      FT_D = 0.0
      IF ((TMAX_D - T) .GT. 0.0 .AND. (T - TMIN_D) .GT. 0.0 .AND.
     &    (TMAX_D - TOT_D) .GT. 0.0 .AND. (TOT_D - TMIN_D) .GT. 0.0)
     &    THEN
          FT_D = ((TMAX_D-T)/(TMAX_D-TOT_D)) *
     &    ((T-TMIN_D)/(TOT_D-TMIN_D))**((TOT_D-TMIN_D)/(TMAX_D-TOT_D))
      END IF
 
      FT_G = MIN(MAX(FT_G, 0.0), 1.0)
      FT_D = MIN(MAX(FT_D, 0.0), 1.0)
      FT   = MIN(MAX(FT_G * FT_D, 0.0), 1.0)
 
      RETURN
      END SUBROUTINE T_DEV
 
!----- Infection rate --------------------------------------------------
      SUBROUTINE F_IR(FT, LWD, Ymax, A, B, IR)
      IMPLICIT NONE
      REAL Ymax, A, B, IR, LWD, FT, ARG
 
      ARG = MAX(A * LWD, 0.0)
      IF (ARG .LE. 0.0) THEN
          IR = 0.0
      ELSE
          IR = Ymax * FT * (1.0 - EXP(-ARG**MAX(B, 1.0E-6)))
      END IF
      IR = MAX(IR, 0.0)
 
      RETURN
      END SUBROUTINE F_IR
 
!----- Fraction of the airborne cloud deposited today ------------------
!  Two ceilings, whichever bites first:
!    settling : only DEP_FRAC of what is airborne reaches the canopy
!               in a day, so a cloud drains over several days;
!    capacity : the canopy holds at most HEALTH_LAI / LESION_S
!               lesion-sized deposition sites.
!  So F_CANSPO x CLOUD = MIN(DEP_FRAC x CLOUD, capacity).
      SUBROUTINE F_CANSPO(LAI, CLOUD, LESION_S, DEP_FRAC, FSS)
      IMPLICIT NONE
      REAL FSS, LAI, CLOUD, LESION_S, DEP_FRAC, CAPACITY
 
      IF (CLOUD .LE. 0.0 .OR. LESION_S .LE. 0.0) THEN
          FSS = 0.0
      ELSE
          CAPACITY = MAX(LAI, 0.0) / LESION_S
          FSS = MIN(CAPACITY / CLOUD, MAX(DEP_FRAC, 0.0))
          FSS = MAX(FSS, 0.0)
      END IF
 
      RETURN
      END SUBROUTINE F_CANSPO
 
!----- Deposited spores ------------------------------------------------
      SUBROUTINE F_DS(FSS, CLOUD, DS)
      IMPLICIT NONE
      REAL FSS, CLOUD, DS
      DS = MAX(FSS * CLOUD, 0.0)
      RETURN
      END SUBROUTINE F_DS
 
!----- Latency progress rate -------------------------------------------
      SUBROUTINE F_LR(FT_D, LDmin, LR)
      IMPLICIT NONE
      REAL LR, FT_D, LDmin
      LR = MAX(FT_D, 0.0) / MAX(LDmin, 1.0E-6)
      RETURN
      END SUBROUTINE F_LR
 
!----- New latent lesions ----------------------------------------------
      SUBROUTINE F_LS(IR, DS, LS)
      IMPLICIT NONE
      REAL IR, DS, LS
      LS = MAX(IR * DS, 0.0)
      RETURN
      END SUBROUTINE F_LS
 
!----- Lesion ageing rate ----------------------------------------------
      SUBROUTINE F_LAR(FT_D, LESLIFEMAX, Lesion_Rate)
      IMPLICIT NONE
      REAL FT_D, LESLIFEMAX, Lesion_Rate
 
      IF (FT_D .LE. 0.0) THEN
          Lesion_Rate = 0.0
      ELSE
          Lesion_Rate = FT_D / MAX(LESLIFEMAX, 1.0E-6)
      END IF
 
      RETURN
      END SUBROUTINE F_LAR
 
!----- Lesion age factor (triangular, peak at LESIONAGEOPT) ------------
!  LA and LESIONAGEOPT are both RELATIVE ages in (0,1), i.e. fractions
!  of LESLIFEMAX. 
      SUBROUTINE F_LAF(LA, LAF, LESIONAGEOPT)
      IMPLICIT NONE
      REAL LAF, LA, LESIONAGEOPT, TMP
 
      IF (LA .LT. LESIONAGEOPT) THEN
          TMP = LA / MAX(LESIONAGEOPT, 1.0E-6)
      ELSE
          TMP = 1.0 - (LA - LESIONAGEOPT)
     &          / MAX(1.0 - LESIONAGEOPT, 1.0E-6)
      END IF
      LAF = MIN(MAX(TMP, 0.0), 1.0)
 
      RETURN
      END SUBROUTINE F_LAF
 
!----- Expanded fraction of a lesion's final area ----------------------
!  Rising limb of the same triangular lesion-age curve F_LAF uses: a
!  lesion breaks latency as a fleck and reaches its final size at
!  LESIONAGEOPT, the age at which the model already places peak
!  sporulation. 
      SUBROUTINE F_LEXP(LA, LESIONAGEOPT, EXPN)
      IMPLICIT NONE
      REAL LA, LESIONAGEOPT, EXPN

      EXPN = MAX(LA, 0.0) / MAX(LESIONAGEOPT, 1.0E-6)
      EXPN = MIN(MAX(EXPN, 0.0), 1.0)

      RETURN
      END SUBROUTINE F_LEXP

!----- Logistic sporulation in population space ------------------------

      SUBROUTINE F_PPSR_POP(KVERHULST, RVERHULST, NPREV, INF_SURF_PREV,
     &                      LAI_SUSC, POT_SPO_PER_AREA)
      IMPLICIT NONE
      REAL KVERHULST, RVERHULST, NPREV, INF_SURF_PREV, LAI_SUSC
      REAL POT_SPO_PER_AREA, KTOTAL, DN, NEW_SPORES
      REAL, PARAMETER :: EPSL = 1.0E-6
 
      POT_SPO_PER_AREA = 0.0
      IF (LAI_SUSC .LE. 0.0) RETURN
      IF (NPREV .LE. 0.0) RETURN
      IF (INF_SURF_PREV .LE. 0.0) RETURN
 
      KTOTAL     = MAX(KVERHULST * LAI_SUSC, EPSL)
      DN         = RVERHULST * NPREV * (KTOTAL - NPREV) / KTOTAL
      NEW_SPORES = MAX(DN, 0.0)
      POT_SPO_PER_AREA = NEW_SPORES / INF_SURF_PREV
 
      RETURN
      END SUBROUTINE F_PPSR_POP
 
!----- DVIP daily infection-probability class (Beruski et al. 2020) ----
      SUBROUTINE CALC_DVIP(LWD, T, DVIP)
      IMPLICIT NONE
      REAL,    INTENT(IN)  :: LWD, T
      INTEGER, INTENT(OUT) :: DVIP
 
      REAL NLes, TERM_LWD, TERM_TEMP, EXP_VAL
      REAL, PARAMETER :: A_COEFF = 12.611
      REAL, PARAMETER :: LWD_OPT = 20.0
      REAL, PARAMETER :: LWD_SIG =  9.0
      REAL, PARAMETER :: T_OPT   = 23.0
      REAL, PARAMETER :: T_SIG   =  5.0
 
      DVIP = 0
      NLes = 0.0
      IF (LWD .LE. 1.0) RETURN
 
      TERM_LWD  = ((LWD - LWD_OPT) / LWD_SIG)**2.0
      TERM_TEMP = ((T   - T_OPT)   / T_SIG)**2.0
      EXP_VAL   = -2.5 * (TERM_LWD + TERM_TEMP)
 
      IF (EXP_VAL .LT. -20.0) THEN
          NLes = 0.0
      ELSE
          NLes = A_COEFF * EXP(EXP_VAL)
      END IF
 
      IF (NLes .LE. 0.5) THEN
          DVIP = 0
      ELSEIF (NLes .LE. 3.0) THEN
          DVIP = 1
      ELSEIF (NLes .LE. 6.0) THEN
          DVIP = 2
      ELSE
          DVIP = 3
      END IF
 
      RETURN
      END SUBROUTINE CALC_DVIP
 
!----- Fungicide decision / application --------------------------------
!  NSprays is the number of applications the model estimates from the
!  risk index.  It is independent of USE_FUNGICIDE, which only decides
!  whether the applications reduce the infection rate.  This lets a
!  user obtain the estimated spray schedule without imposing its
!  effect on the epidemic.
      SUBROUTINE APPLY_FUNGICIDE(DVIP_today, DVIP_pts, idx, SUM7,
     &    BufferDays, FungActive, ResidualDays, NSprays,
     &    USE_FUNGICIDE, FUNG_RES_D, FUNG_BUF_D, DVIP_THR, HEALTH_LAI)
 
      IMPLICIT NONE
      INTEGER DVIP_today, DVIP_pts(7), idx, SUM7
      INTEGER BufferDays, ResidualDays, NSprays
      INTEGER FUNG_RES_D, FUNG_BUF_D, DVIP_THR
      LOGICAL FungActive, USE_FUNGICIDE
      REAL    HEALTH_LAI
      LOGICAL CAN_SPRAY
 
!     1. Rolling 7-day risk sum
      SUM7 = SUM7 - DVIP_pts(idx) + DVIP_today
      DVIP_pts(idx) = DVIP_today
      idx = MOD(idx, 7) + 1
 
!     2. Re-spray buffer
      IF (BufferDays .GT. 0) BufferDays = BufferDays - 1
 
!     3. Residual protection
      IF (FungActive) THEN
          ResidualDays = ResidualDays - 1
          IF (ResidualDays .LE. 0) FungActive = .FALSE.
      END IF
 
!     4. Decision.  
      CAN_SPRAY = (HEALTH_LAI .GE. 0.1)
      IF (SUM7 .GE. DVIP_THR .AND. BufferDays .EQ. 0 .AND. CAN_SPRAY)
     &    THEN
          NSprays    = NSprays + 1
          BufferDays = MAX(FUNG_BUF_D, 1)
          IF (USE_FUNGICIDE) THEN
              FungActive   = .TRUE.
              ResidualDays = MAX(FUNG_RES_D, 1)
          END IF
      END IF
 
      RETURN
      END SUBROUTINE APPLY_FUNGICIDE
 
!----- Virtual lesions (Bastiaans 1991; Primiano & Amorim 2020) --------

      SUBROUTINE F_VIRTUAL_LESIONS(s, beta, fvl)
      IMPLICIT NONE
      REAL s, beta, fvl, SS, BSAFE
 
      SS    = MIN(MAX(s, 0.0), 0.999999)
      BSAFE = MAX(beta, 1.0)
 
      IF (BSAFE .LE. 1.000001) THEN
          fvl = 1.0
      ELSE
          fvl = (1.0 - SS) ** (BSAFE - 1.0)
      END IF
      fvl = MIN(MAX(fvl, 0.0), 1.0)
 
      RETURN
      END SUBROUTINE F_VIRTUAL_LESIONS
 
!----- Severity --------------------------------------------------------

      SUBROUTINE F_SEVERITY(LAI_NECRO, LAI_PEAK, SEV)
      IMPLICIT NONE
      REAL, INTENT(IN)  :: LAI_NECRO, LAI_PEAK
      REAL, INTENT(OUT) :: SEV
 
      SEV = 100.0 * MAX(LAI_NECRO, 0.0) / MAX(LAI_PEAK, 1.0E-6)
      SEV = MIN(MAX(SEV, 0.0), 100.0)
 
      RETURN
      END SUBROUTINE F_SEVERITY
 
!----- Disease-induced senescence --------------------------------------
      SUBROUTINE F_DEFOLIATION(WTLF, SLDOT, SEVERITY_PCT, rrds,
     &                         DISEASE_SEN_RATE)
      IMPLICIT NONE
      REAL, INTENT(IN)  :: WTLF, SLDOT, SEVERITY_PCT, rrds
      REAL, INTENT(OUT) :: DISEASE_SEN_RATE
 
      REAL rrsen, rrsenD
      REAL, PARAMETER :: EPS = 1.0E-6
 
      IF (WTLF .GT. EPS) THEN
!         Relative senescence rate already imposed by DSSAT
          rrsen  = SLDOT / WTLF
!         Relative senescence rate attributable to disease
          rrsenD = rrds * (MIN(MAX(SEVERITY_PCT,0.0),100.0) / 100.0)
!         Mass to remove, avoiding double counting on the same area
          DISEASE_SEN_RATE = (rrsenD - (rrsen * rrsenD)) * WTLF
          DISEASE_SEN_RATE = MAX(0.0, DISEASE_SEN_RATE)
      ELSE
          DISEASE_SEN_RATE = 0.0
      END IF
 
      RETURN
      END SUBROUTINE F_DEFOLIATION
 
!
!***********************************************************************
!  VARIABLE LISTING
!***********************************************************************
! --------------------------- Arguments --------------------------------
! DYNAMIC           : DSSAT phase flag (RUNINIT/SEASINIT/RATE/OUTPUT/
!                     SEASEND)
! CONTROL           : ControlType (%DAS, %RUN, %YRSIM, %FILEIO, ...)
! ISWITCH           : SwitchType (not used; kept for interface parity)
! Tmin, Tmax (C)    : Daily min/max air temperature
! RH (%)            : Relative humidity
! LAI_TOTAL (m2 m-2): Canopy leaf area index from the crop model
! WTLF (kg ha-1)    : Leaf mass
! SLDOT (kg ha-1 d-1): Natural leaf senescence rate from DSSAT
! ESP_LAT_HIST(:,:) : Cohort ring buffer (MAXDAYS x 5); columns
!                     C_LES / C_LATP / C_INFF / C_AGE
! YRDOY, YREMRG     : Current date, emergence date (YRDOY convention)
! NVEG0             : DAS threshold for the emergence gate
! YREND             : Harvest / end-of-season date
! DISEASE_LAI       : OUTPUT, cumulative diseased area (cm2 m-2)
! VIRTUAL_PHOTO_FACTOR : OUTPUT (0..1), multiplies PGAVL in CROPGRO
! DISEASE_SEN_RATE  : OUTPUT (kg ha-1 d-1), disease-induced senescence
!
! --------------------------- Fixed limits -----------------------------
! MAXDAYS      : Cohort ring buffer length (days)
! MAXPRESEASON : Maximum pre-plant days stored for output replay
! EPS          : Small epsilon for safe divisions
!
! ------------- Parameters (parameter file, 25 tokens) -----------------
! Token numbers are the record columns READ_DISEASE_PARAMETERS reads.
! Tokens 1-2 are VAR# and VRNAME.  There are no optional columns: a
! record with fewer than 25 tokens is rejected.
!
!  3 DAE_START        : Observed onset date, in days after emergence,
!                     when the user has one from the field.  Arrival
!                     then happens on that day, calendarised, and
!                     FAV_THR is not read.  -99 = no observed date,
!                     onset is left to the weather.  Held internally as
!                     DAE_MIN / DAE_MIN_PRESENT.
!  4 NCYCLE           : 'P' polycyclic, 'M' monocyclic
!  5 SRC_HALF         : Half-saturation of the source response, on the
!                     SOURCE_PRESSURE scale
!  6 FAV_THR          : Accumulated favourability at which external
!                     inoculum arrives, in favourable-day equivalents
!                     on the FAV_SUM scale.  ONSET TIMING when
!                     DAE_START is -99; ignored otherwise.
!  7 NDS              : Number of depositable primary spores released on
!                     arrival, before the F_SOURCE seasonal modifier.
!                     ONSET MAGNITUDE.  Cannot move the arrival day.
!                     Same role it always had -- the primary inoculum
!                     scale -- but now the whole dose delivered once,
!                     not a per-day rate against the cloud index, so it
!                     needs recalibrating on the new scale.
!  8 DEP_FRAC         : Fraction of the airborne cloud that settles
!                     onto the canopy in one day, in (0,1].  Caps
!                     F_CANSPO, so a cloud drains over several days
!                     instead of all at once; 1.0 restores the old
!                     instantaneous deposition.  Applies to BOTH the
!                     primary pool and the secondary cloud.
!                     NOT a free coordinate for the calibration: it
!                     is partly collinear with NDS on the primary
!                     path and with RVERHULST on the secondary one.
!                     Fix it from literature or a predeclared
!                     sensitivity analysis, the way SRC_HALF is
!                     fixed by diagnostic.
!  9 LESION_S (m2)    : Average lesion surface area
! 10 KVERHULST        : Logistic carrying capacity per unit LAI
! 11 RVERHULST        : Logistic intrinsic rate
! 12 YMAX             : Maximum infection efficiency (0-1)
! 13 COF_A, 14 COF_B  : Infection response vs leaf wetness
! 15 TMIN_G 16 TOT_G 17 TMAX_G : Cardinal temperatures, germination
! 18 TMIN_D 19 TOT_D 20 TMAX_D : Cardinal temperatures, development
! 21 LDMIN (d)        : Minimum latent period
! 22 LESIONAGEOPT     : Relative lesion age of peak sporulation, in (0,1)
! 23 LESLIFEMAX (d)   : Maximum lesion lifespan
! 24 BETA             : Virtual lesion exponent; 1 = no effect
! 25 RRDS (d-1)       : Relative disease senescence rate
!
! ------------------ Fixed constants (set in RUNINIT) ------------------
! LAI_MIN_START     : Minimum LAI for deposition
! SPOR_DECAY, SEC_DECAY, SRC_SURV : Daily retention factors of the
!                     primary pool, the secondary cloud and the regional
!                     source (fixed in RUNINIT)
! USE_FUNGICIDE, FUNG_EFFICIENCY, FUNG_RES_D, FUNG_BUF_D, DVIP_THR :
!                     Fungicide block, fixed in RUNINIT; belongs in the
!                     FILEX management section
! INOC_LOOKBACK (d) : Pre-plant replay window (fixed)
!
! ------------------------- Environmental state ------------------------
! DAILY_IP (0..1)   : Daily climatic favourability for the source
! SOURCE_PRESSURE   : Slow decayed sum of DAILY_IP (34-day half-life);
!                     through F_SOURCE it scales the arrival dose
! FAV_SUM           : Undecayed sum of DAILY_IP from EMERGENCE onward;
!                     the arrival clock.  NOT counted during the
!                     pre-plant replay -- see COUNT_EXPOSURE
! COUNT_EXPOSURE    : .TRUE. only once the crop has emerged; gates
!                     FAV_SUM so the clock measures the CROP's exposure
!                     and not the site's pre-season climate
! PRI_RELEASED      : Arrival latch, fires once per season
! PRI_POOL          : Primary inoculum still airborne after arrival;
!                     decays at SPOR_DECAY and is depleted by DS_PRI
! SEC_SPORE_CLOUD   : Airborne secondary inoculum
! SEC_SPORES_PENDING: Today's emission, airborne tomorrow
!
! --------------------------- Epidemic state ---------------------------
! LAI_PEAK_SEASON   : Maximum LAI observed this season
! CUM_NECROTIC      : Cumulative necrotic LAI (m2 m-2)
! LAI_SUSC          : LAI_PEAK_SEASON - CUM_NECROTIC (epidemic ref.)
! HEALTH_LAI        : LAI_TOTAL - CUM_NECROTIC (green tissue present)
! PREV_IS           : Infectious surface yesterday
! SLOT_DAE(:)       : DAE that created each ring slot; 0 = free
! N_ACTIVE_COH      : Number of occupied ring slots
! ADMITTED_AREA(:)  : Area admitted to each cohort at activation
! SEVERITY_PCT      : Severity (%), output column SEV%; drives
!                     defoliation and is the calibration target
! LAI_SUSC          : epidemic reference -- tissue that ever existed minus
!                     what the epidemic already destroyed.  Both terms are on
!                     the PEAK basis, so the subtraction is well posed.  Used
!                     for area admission, carrying capacity and severity.
! HEALTH_LAI        : the part of TODAY'S canopy that is not diseased.  Used
!                     for deposition -- spores only land on what is there.
!
! --------------- Output columns (DISMO.OUT) ---------------
! Written in this order.  Pre-plant rows replay the values the weather
! had during the DISMO_PRESEASON reconstruction, so RHU%, LWDh, FTMP and
! FAVS are meaningful before planting; the epidemic columns are zero.
!
! YEAR DOY : calendar date        DAS  : days since simulation start,
!                                        which is planting - INOC_LOOKBACK
! DAE  : days after emergence, 0 before the crop is up
! LAIH : HEALTH_LAI       LWDh : leaf wetness   RHU% : humidity
! FTMP : temperature response          LAIT : LAI_TOTAL
! SUM7 : 7-day DVIP sum   NSPRAYS : estimated applications
! FACT : fungicide residual active     SEV% : severity, for calibration
! SRCP : SOURCE_PRESSURE
! FAVS : FAV_SUM, the arrival clock -- read it on the observed date of
!        first symptom to estimate FAV_THR directly from the data
! PRIM : PRI_CLOUD, primary inoculum available today; the first non-zero
!        day is the arrival day
! SSCL : SEC_SPORE_CLOUD
! DSTO : deposited spores LSTO : new lesions   NLTO : new necrosis
! FSRC : source response, scales the dose released on arrival
!
! --------------------------- Subroutines ------------------------------
! DISEASE_LEAF       : The module itself, called once per DYNAMIC phase
! DISMO_SEASON_RESET : Clears all season-scope state (single source)
! DISMO_CLEAR_SLOT   : Frees one cohort ring slot
! DISMO_SPLIT        : Whitespace tokeniser for the parameter file
! DISMO_TOKR         : Token to REAL, with default and -99 handling
! READ_DISEASE_PARAMETERS : Parameter file reader
! DISMO_PRESEASON    : Pre-plant weather replay
! DISMO_UPDATE_ENVIRONMENT : Daily T, LWD, temperature responses and
!                      inoculum accumulators
! DISMO_NORMALIZE_DATE / DISMO_SHIFT_DATE / DISMO_READ_WTH_FIELD
! F_TAVG / F_DEW / F_RH / F_LWD : Weather derivations
! T_DEV              : Beta temperature responses (FT_G, FT_D, FT)
! F_IR               : Infection rate vs FT and LWD
! F_CANSPO / F_DS    : Canopy interception capacity and deposition
! F_LR / F_LS        : Latency rate and new latent lesions
! F_LAR / F_LAF      : Lesion ageing rate and lesion age factor
! F_LEXP             : Expanded fraction of a lesion's final area
! F_PPSR_POP         : Logistic sporulation (see note B below)
! CALC_DVIP / APPLY_FUNGICIDE : Risk index and spray decision
! F_VIRTUAL_LESIONS  : Residual photosynthetic reduction
! F_SEVERITY         : Severity (%)
! F_DEFOLIATION      : Disease-induced senescence rate
!                   
! ------------------------- Arrival and Onset logic ---------------------
!  Arrival is an EVENT that opens a TAP, not a single shot.  The latch
!  fires once per season, on the day the weather clock (or the observed
!  onset) says the pathogen reached the field; from that day on the
!  regional source keeps supplying, modulated by today's favourability.
!
!  SOURCE_PRESSURE is already carried all season with its own
!  34-day half-life; this makes the dose it delivers a daily quantity,
!  gated by the same DAILY_IP that builds FAV_SUM, so external supply
!  follows the weather instead of a single calendar day.
!
!  Two ways to decide the day, chosen by DAE_START in the input file:
!
!  DAE_START = a value : the user knows the onset date from the field.
!                        Arrival happens on that day after emergence,
!                        calendarised, and the weather clock is ignored.
!  DAE_START = -99     : no observed date, so onset is decided by the
!                        weather.  Arrival happens when FAV_SUM, the
!                        accumulated favourability counted from the
!                        first replayed pre-plant day, reaches FAV_THR.
!
!  FAV_THR sets WHEN and NDS sets HOW MUCH.  They act on different
!  features of the severity curve -- NDS cannot move the arrival
!  day and FAV_THR cannot change the size of the supply -- so they are
!  separately identifiable against a severity objective.  NDS is the
!  dose delivered on a FULLY FAVOURABLE day (DAILY_IP = 1) rather than
!  the whole season's dose, so values fitted under the old meaning do
!  not carry over and must be re-estimated.
!     
!     SOURCE_PRESSURE : slow store, 34-day half-life, steady state
!                       ~ DAILY_IP / (1 - SRC_SURV).  A decayed SUM: it
!                       measures how large the regional inoculum source
!                       has grown, and through F_SOURCE it scales the
!                       dose released on arrival.  Calibrate SRC_HALF on
!                       THIS scale.
!     FAV_SUM         : undecayed running sum of DAILY_IP, in
!                       favourable-day equivalents, counted FROM
!                       EMERGENCE ONWARD (COUNT_EXPOSURE).  It is the
!                       clock the arrival latch reads, and the only
!                       DISMO state that carries no decay -- a crop is
!                       either far enough into favourable weather for
!                       the pathogen to have arrived, or it is not.
!=======================================================================
