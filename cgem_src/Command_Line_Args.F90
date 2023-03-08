      Subroutine Command_Line_Args(ivar)

      use cgem_vars

      implicit none

      integer, intent(out) :: ivar
      character(10) :: print_var
      integer c_count

! ./CGEM Which_code InputFile InitializationFile OutputFileBase 

! --- Command Line Arguments for file names ---
         c_count = command_argument_count()

       if (c_count.gt.0) then
       call get_command_argument(1,print_var)  !User selects what to print
       select case (print_var)
       case ("A")
         ivar = iA(1)
       case ("Qn")
         ivar = iQn(1)
       case ("Qp")
         ivar = iQp(1)
       case ("Z")
         ivar = iZ(1)
       case ("NO3")
         ivar = iNO3
       case ("NH4")
         ivar = iNH4
       case ("PO4")
         ivar = iPO4
       case ("DIC")
         ivar = iDIC
       case ("O2")
         ivar = iO2
       case ("OM1_A")
         ivar = iOM1_A
       case ("OM2_A")
         ivar = iOM2_A
       case ("OM1_Z")
         ivar = iOM1_Z
       case ("OM2_Z")
         ivar = iOM2_Z
       case ("OM1_R")
         ivar = iOM1_Z
       case ("OM2_R")
         ivar = iOM2_Z
       case ("OM1_BC")
         ivar = iOM1_BC
       case ("OM2_BC")
         ivar = iOM2_BC
       case ("CDOM")
         ivar = iCDOM
       case ("Si")
         ivar = iSi
       case ("Alk")
         ivar = iAlk
       case ("Tr")
         ivar = iTr

       case default
        write(6,*) "These are the only options, case sensitive:"
        write(6,*) "A,Qn,Qp,Z,NO3,NH4,PO4,DIC,O2"
        write(6,*) "OM1_A,OM2_A,OM1_Z,OM2_Z,OM1_R,OM2_R"
        write(6,*) "CDOM,Si,Alk,Tr"
        write(6,*) "or you can choose nothing."
        write(6,*) "You chose ", print_var
        STOP
      end select

      else
        ivar=0
      endif

!-A; Phytoplankton number density (cells/m3);
!-Qn: Phytoplankton Nitrogen Quota (mmol-N/cell)
!-Qp: Phytoplankton Phosphorus Quota (mmol-P/cell)
!-Z: Zooplankton number density (individuals/m3);
!-NO3; Nitrate (mmol-N/m3)
!-NH4; Ammonium (mmol-N/m3)
!-PO4: Phosphate (mmol-P/m3)
!-DIC: Dissolved Inorganic Carbon (mmol-C/m3) 
!-O2: Molecular Oxygen (mmol-O2/m3)
!-OM1_A: (mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           dead Phytoplankton
!-OM2_A: (mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!           dead Phytoplankton 
!-OM1_Z:(mmol-C/m3--particulate)
!        -- Particulate Organic Matter arising from 
!           Zooplankton fecal pellets.
!-OM2_Z:(mmol-C/m3--dissolved)
!        -- Dissolved Organic Matter arising from 
!          Zooplankton fecal pellets.
!-OM1_R: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter arising from river outflow
!-OM2_R: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter arising from river outflow
!-CDOM: (ppb) 
!        -- Colored Dissolved Organic Matter
!-Silica: (mmol-Si/m3) 
!        -- Silica
!-OM1_BC: (mmol-C/m3--particulate)
!         -- Particulate Organic Matter in initial and boundary 
!            conditions 
!-OM2_BC: (mmol-C/m3--dissolved)
!         -- Dissolved Organic Matter in initial and boundary
!            conditions
!-ALK:  (mmol-HCO3/m3)?
!        -- Alkalinity
!- Tr: Tracer
 

       !if (c_count.gt.1) then
       !!  call get_command_argument(2,input_filename)  !User selects input file name
       !endif

       !if (c_count.gt.2) then
       !!  call get_command_argument(3,init_filename) !User selects initial conditions file name
       !endif

       return

       END Subroutine Command_Line_Args
