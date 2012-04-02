MODULE geometry

  USE Trianlge_Class
  USE Quadrangle_Class
  
  IMPLICIT NONE
  PRIVATE

  !==========================================================
  TYPE :: type_verts
     INTEGER, DIMENSION(:), POINTER :: verts
     INTEGER, DIMENSION(:), POINTER :: NU_seg
     INTEGER, DIMENSION(:), POINTER :: n_ele
  END TYPE type_verts
  TYPE(type_verts), DIMENSION(:), ALLOCATABLE :: ele
  !----------------------------------------------------------
  
  TYPE :: elements_ptr
     CLASS(element), POINTER :: p
  END type elements_ptr
  TYPE(elements_ptr), DIMENSION(:), ALLOCATABLE :: elements

  !----------------------------------------------------------

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: edge_ele
  !==========================================================
  
  INTEGER :: N_nodes, N_ele_mesh

  REAL(KINd=8), DIMENSION(:,:), ALLOCATABLE :: rr_nodes

  INTEGER, PARAMETER :: N_dim = 2

  INTEGER :: N_dofs, N_elements, N_seg
  !----------------------------------------------------------
  
  ! Gmsh standard element types
  INTEGER, PARAMETER :: GMSH_LINE       = 1
  INTEGER, PARAMETER :: GMSH_TRIANGLE   = 2
  INTEGER, PARAMETER :: GMSH_QUADRANGLE = 3
  INTEGER, PARAMETER :: GMSH_POINT      = 15
  !==========================================================

  PUBLIC :: read_gmshMesh, read_MeshFBx, Init_Elements
  PUBLIC :: N_dim, N_dofs, N_elements, elements
  !==========================================================

CONTAINS

  !==============================
  SUBROUTINE Init_Elements(Order)
  !==============================

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Order
    !---------------------------------------------------

    TYPE(triangle),   POINTER :: tri
    TYPE(quadrangle), POINTER :: qua

    INTEGER,      DIMENSION(:),    ALLOCATABLE :: VV
    REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: RR

    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_2_q 

    INTEGER :: Nv, i, j, jq, jt, is1, is2, istat
    !-----------------------------------------------------

    jq = 0

    CALL find_segments()

    SELECT CASE(Order)

    !==========================================================
    CASE(2)
    !==========================================================

       N_dofs     = N_nodes       
       N_elements = N_ele_mesh

       ALLOCATE( elements(N_elements) )
       
       DO jt = 1, N_elements

          Nv = SIZE(ele(jt)%verts(:))

          SELECT CASE(Nv)

          !--------------------
          CASE(3) ! TRIANGLE P1
          !----------------------------------------------------

             ALLOCATE(tri, STAT=istat)
             IF(istat /=0) THEN
                WRITE(*,*) 'ERROR: failed trianlge allocation'
             ENDIF
             CALL tri%initialize( "element", ele(jt)%verts,     &
                                  rr_nodes(:, ele(jt)%verts),   &
                                  ele(jt)%NU_seg, ele(jt)%n_ele )

             elements(jt)%p => tri

          !----------------------
          CASE(4) ! QUADRANGLE Q1
          !----------------------------------------------------
             
             ALLOCATE(qua, STAT=istat)
             IF(istat /=0) THEN
                WRITE(*,*) 'ERROR: failed trianlge allocation'
             ENDIF
             CALL qua%initialize( "element", ele(jt)%verts,     &
                                  rr_nodes(:, ele(jt)%verts),   &
                                  ele(jt)%NU_seg, ele(jt)%n_ele )

             elements(jt)%p => qua

          END SELECT
          
       ENDDO

    !==========================================================
    CASE(3)
    !==========================================================

       N_dofs     = N_nodes + N_seg       
       N_elements = N_ele_mesh

       ALLOCATE( elements(N_elements) )
       
       DO jt = 1, N_elements

          Nv = SIZE(ele(jt)%verts(:))

          SELECT CASE(Nv)

          !--------------------
          CASE(3) ! TRIANGLE P2
          !----------------------------------------------------

             ALLOCATE( VV(6), RR(N_dim, 6) )

             VV(1:3) = ele(jt)%verts(:)
             VV(4:6) = ele(jt)%NU_seg( (/3,1,2/)  ) + N_nodes

             DO i = 1, Nv

                RR(:, i) = rr_nodes(:, ele(jt)%verts(i))

                j = MOD(i, MIN(i+1, Nv)) + 1

                is1 = ele(jt)%verts(i) 
                is2 = ele(jt)%verts(j) 

                RR(:, Nv+i) = 0.5d0*( rr_nodes(:, is1) + rr_nodes(:, is2) )

             ENDDO

             ALLOCATE(tri, STAT=istat)
             IF(istat /=0) THEN
                WRITE(*,*) 'ERROR: failed trianlge allocation'
             ENDIF
             CALL tri%initialize( "element", VV, RR,            &
                                  ele(jt)%NU_seg, ele(jt)%n_ele )

             elements(jt)%p => tri

             DEALLOCATE( VV, RR )

          !----------------------
          CASE(4) ! QUADRANGLE Q2
          !----------------------------------------------------

             jq = jq + 1

             ALLOCATE( VV(9), RR(N_dim, 9) )

             VV(1:4) = ele(jt)%verts(:)
             VV(5:8) = ele(jt)%NU_seg( (/4,1,2,3/) ) + N_nodes
             VV(9)   = jq + N_nodes + N_seg

             DO i = 1, Nv

                RR(:, i) = rr_nodes(:, ele(jt)%verts(i))

                j = MOD(i, MIN(i+1, Nv)) + 1

                is1 = ele(jt)%verts(i)
                is2 = ele(jt)%verts(j)
               
                RR(:, Nv+i) = 0.5d0 * ( rr_nodes(:, is1) + rr_nodes(:, is2) )

             ENDDO

             RR(1, 9) = SUM( RR(1, 1:Nv) ) / REAL(Nv)
             RR(2, 9) = SUM( RR(2, 1:Nv) ) / REAL(Nv)                
             
             ALLOCATE(qua, STAT=istat)
             IF(istat /=0) THEN
                WRITE(*,*) 'ERROR: failed trianlge allocation'
             ENDIF
             CALL qua%initialize( "element", VV, RR,             &
                                  ele(jt)%NU_seg, ele(jt)%n_ele ) 

             elements(jt)%p => qua

             DEALLOCATE( VV, RR )

          END SELECT
          
       ENDDO       

    END SELECT

    N_dofs = N_dofs + jq

    !*************************
    CALL loc2loc_connectivity()

    DEALLOCATE( ele, edge_ele, rr_nodes )
    
  END SUBROUTINE Init_Elements
  !===========================

  !================================
  SUBROUTINE loc2loc_connectivity()
  !================================
  !
  ! Attach to each face of the elment the data from 
  ! neighbour elements which share the same face

    IMPLICIT NONE

    REAL(KIND=8), DIMENSION(:,:,:), POINTER :: p_Dphi_2_q

    INTEGER :: jt, j, n_ele, i, iq, jq, k1, k2
    !-----------------------------------------


    DO jt = 1, N_elements

       DO j = 1, elements(jt)%p%N_faces

          n_ele = elements(jt)%p%faces(j)%f%c_ele

          IF( n_ele /= 0) THEN

             ! Find the local numeration of the face
             DO i = 1, elements(n_ele)%p%N_faces
                
                IF( elements(n_ele)%p%faces(i)%f%g_seg == &
                       elements(jt)%p%faces(j)%f%g_seg ) EXIT
                
             ENDDO

             p_Dphi_2_q => elements(n_ele)%p%faces(i)%f%p_Dphi_1_q

             ! Gradient of basis functions
             ALLOCATE( elements(jt)%p%faces(j)%f%p_Dphi_2_q( &
                         SIZE(p_Dphi_2_q, 1), &
                         SIZE(p_Dphi_2_q, 2), &
                         SIZE(p_Dphi_2_q, 3)  ) &
                     )

             ! Gradient of the solution        
             ALLOCATE( elements(jt)%p%faces(j)%f%p_Du_2_q( &
                         SIZE(p_Dphi_2_q, 1), &
                         SIZE(p_Dphi_2_q, 3)  ) &
                     )
                                            
             DO jq = 1, elements(jt)%p%faces(j)%f%N_quad

                ! Find the local numeration of the quadrature
                ! points
                DO iq = 1, elements(n_ele)%p%faces(i)%f%N_quad

                   IF( SUM( elements(n_ele)%p%faces(i)%f%xx_q(:, iq) - &
                               elements(jt)%p%faces(j)%f%xx_q(:, jq) ) == 0.d0 ) EXIT

                ENDDO
                             
                elements(jt)%p%faces(j)%f%p_Dphi_2_q(:,:, jq) = p_Dphi_2_q(:,:, iq)

             ENDDO             

             NULLIFY( p_Dphi_2_q )

             !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
             ! Correspondence between the local numeration of the DOFs\
             ! of two neighbouring elements which share the face 'j'  \
             !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
             ALLOCATE( elements(jt)%p%faces(j)%f%loc_con( &
                       elements(jt)%p%N_points) &
                     )

             elements(jt)%p%faces(j)%f%loc_con = 0

             DO k1 = 1, elements(jt)%p%N_points

                DO k2 = 1, elements(n_ele)%p%N_points

                   IF( SUM( elements(n_ele)%p%Coords(:, k2) - &
                               elements(jt)%p%Coords(:, k1) ) == 0.d0) THEN

                      elements(jt)%p%faces(j)%f%loc_con(k1) = k2
                      
                   ENDIF

                ENDDO
                
             ENDDO
             !\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
           
             
          ELSE ! No neighbouring element (boundary face)

             p_Dphi_2_q => elements(jt)%p%faces(j)%f%p_Dphi_1_q

             ! Gradient of basis functions             
             ALLOCATE( elements(jt)%p%faces(j)%f%p_Dphi_2_q( &
                         SIZE(p_Dphi_2_q, 1), &
                         SIZE(p_Dphi_2_q, 2), &
                         SIZE(p_Dphi_2_q, 3)  ) &
                      )

             ! Gradient of the solution        
             ALLOCATE( elements(jt)%p%faces(j)%f%p_Du_2_q( &
                         SIZE(p_Dphi_2_q, 1), &
                         SIZE(p_Dphi_2_q, 3)  ) &
                     )

             elements(jt)%p%faces(j)%f%p_Dphi_2_q = 0.d0

             ALLOCATE( elements(jt)%p%faces(j)%f%loc_con(elements(jt)%p%N_points) )

             elements(jt)%p%faces(j)%f%loc_con = 0
             
             NULLIFY( p_Dphi_2_q )

          ENDIF
          

       ENDDO       

    ENDDO

  END SUBROUTINE loc2loc_connectivity
  !==================================
  
  !=========================
  SUBROUTINE find_segments()
  !=========================
  !
  ! Algorithm to find the segment on the mensh
  !   

    IMPLICIT NONE

    INTEGER, PARAMETER :: v_max = 15
         
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: J_CON
    INTEGER, DIMENSION(:,:),   ALLOCATABLE :: cn, cn_
    INTEGER, DIMENSION(:,:),   ALLOCATABLE :: edge_ele_
    INTEGER, DIMENSION(:),     ALLOCATABLE :: seg_b

    LOGICAL :: boundary
    
    INTEGER :: N_segb, found_seg, found_segb, ig0
    INTEGER :: is_1, is_2, is_3, is_4, cn_1, cn_2, ele_1, ele_2
    INTEGER :: l_1, l_2, l_3, l_4, Nv, Ns
    INTEGER :: jt, k, j, kv, iseg, i

    !========================================================

    N_seg = N_nodes + N_ele_mesh - 1

    ALLOCATE ( J_CON(N_nodes, v_max, 2) )
          
    ALLOCATE ( cn_(2, N_seg), edge_ele_(2, N_seg) )

    J_CON = 0; cn_ = 0
      
    found_seg = 0; ig0 = 0
      
    edge_ele_ = 0

    !
    !------------------------------------------------------------------------
    !       Segments in the whole domain
    !------------------------------------------------------------------------
    !
    !  edge_ele(:, i) -- Connectivity between the edge i and elements
    !        cn (:,i) -- Couple of nodes of the segment i
    !
    DO jt = 1, N_ele_mesh

       Nv = SIZE(ele(jt)%verts(:))

       DO k = 1, Nv

          i = MOD(k, MIN(k+1, Nv)) + 1
          j = MOD(i, MIN(k+2, Nv)) + 1
          
          is_1 = ele(jt)%verts(i)          
          is_2 = ele(jt)%verts(j)

          DO kv = 1, v_max
                         
             IF ( J_CON(is_1, kv, 1) == 0 ) GOTO 110
                            
              IF ( J_CON(is_1, kv, 1) == is_2) THEN
               
                 ig0 = J_CON(is_1, kv, 2)
                                   
                 cn_1 = ABS(cn_(1, ig0))
                  
                 cn_(1, ig0) = -cn_(1, ig0)
                  
                 cn_2 = cn_(2, ig0)
                  
                 IF ( is_1 == cn_1 .AND. is_2 == cn_2 ) THEN

                    ! Nothing happens here

                 ELSE 
                  
                    IF (is_1 == cn_2 .AND. is_2 == cn_1) THEN

                       ! Nothing happens here
                        
                    ELSE
                                            
                       WRITE(*,*) 'ERROR: incosistent nodes for the segment', is_1, is_2, cn_1, cn_2
                       WRITE(*,*) 'STOP!'
                        
                       STOP
                        
                    ENDIF
                                      
                 ENDIF
                  
                 edge_ele_(2, ig0) =  jt
                          
                 GOTO 210
          
              ENDIF
               
           ENDDO
           
            
110        DO kv = 1, v_max
                          
              IF ( J_CON(is_2, kv, 1) == 0 ) GOTO 130
                             
              IF ( J_CON(is_2, kv, 1) == is_1) THEN
               
                 ig0 = J_CON(is_2, kv, 2)
                                   
                 cn_1 = ABS(cn_(1, ig0))
                  
                 cn_(1, ig0) = -cn_(1, ig0)
                  
                 cn_2 = cn_(2, ig0)
                  
                 IF ( is_1 == cn_1 .AND. is_2 == cn_2 ) THEN
                    
                    ! Nothing happens here
                    
                 ELSE
                                      
                    IF (is_1 == cn_2 .AND. is_2 == cn_1) THEN
                       
                       ! Nothing happens here
                       
                    ELSE
                                            
                       WRITE(*,*) 'ERROR: incosistent nodes for the segment', is_1, is_2, cn_1, cn_2
                       WRITE(*,*) 'STOP!'
                                               
                       STOP           
                        
                    ENDIF                    
                  
                 ENDIF                 
                  
                 edge_ele_(2, ig0) =  jt                 
                  
                 GOTO 210                
               
              ENDIF              
            
           ENDDO
           
            
           kv = kv + 1           
            
           IF (kv > v_max) THEN
                       
              WRITE(*,*) 'ERROR: v_max exceeded'
              WRITE(*,*) 'STOP!'
                          
              STOP

           ENDIF           
            
130        found_seg = found_seg + 1
                    
           IF (found_seg > N_seg) THEN
                       
              WRITE(*,*) 'ERROR: number of segments exceeded'
              WRITE(*,*) 'STOP!'
            
              STOP
                       
           ENDIF
           
           J_CON(is_2, kv, 1) = is_1
           J_CON(is_2, kv, 2) = found_seg
         
           cn_(1, found_seg) = -is_1
           cn_(2, found_seg) =  is_2
         
           ig0 = found_seg
                  
           edge_ele_(1, ig0) =  jt
                    
210        CONTINUE   
         
        ENDDO
         
     ENDDO  
               
     DEALLOCATE (J_CON)

     !-----------------------
     ! True segment structure
     !-----------------------------------------
     N_seg = COUNT(cn_(1, :) /= 0)

     ALLOCATE(       cn(2, N_seg) )
     ALLOCATE( edge_ele(2, N_seg) )

     DO iseg = 1, N_seg

        IF( cn_(1, iseg) /= 0 ) THEN

                 cn(:, iseg) = cn_(:, iseg)
           edge_ele(:, iseg) = edge_ele_(:, iseg)

        ENDIF

     ENDDO
     

     DEALLOCATE (cn_, edge_ele_)

     !
     !---------------------------------------------------------------------------
     !       Segments belongig to the boundaries
     !---------------------------------------------------------------------------
     !
     ! seg_b(j) -- Connectivity between the boundary segment j and the domain segment
     !
     found_segb = 0
      
     DO iseg = 1, N_seg
        IF (cn(1, iseg) < 0) THEN
            
           found_segb = found_segb + 1
                       
        ENDIF        
      
     ENDDO      
      
     ALLOCATE( seg_b(found_segb) )
      
     found_segb = 0
      
     DO iseg = 1, N_seg
      
        IF (cn(1, iseg) < 0) THEN
                    
           found_segb = found_segb + 1
            
           seg_b(found_segb) = iseg
            
           cn(1, iseg) = -cn(1, iseg)
         
        ENDIF
              
     ENDDO
           
     N_segb = found_segb
      
     IF (MINVAL(cn) <= 0) THEN
      
        WRITE(*,*) 'ERROR: Unable to find all bondary edges'
        WRITE(*,*) 'STOP!'
         
        STOP
      
     ENDIF
    
     WRITE(*,*)
     WRITE(*,*) 'Number of domain segments: ',found_seg
     WRITE(*,*) 'Number of boundary segments: ',found_segb

     !
     !---------------------------------------------------------------------------
     !       Connectivity of each node with one edge
     !---------------------------------------------------------------------------
     !
     ! NU_seg(i, jt) -- local-to-glogal segments connectivity
     !
     DO jt = 1, N_ele_mesh
        Ns = SIZE(ele(jt)%verts(:))
        ALLOCATE(ele(jt)%NU_seg(Ns))
     ENDDO     

     DO iseg = 1, N_seg

        ele_1 = edge_ele(1, iseg)
        ele_2 = edge_ele(2, iseg)

        boundary = .FALSE.
         
        IF (ele_1 == 0 .OR. ele_2 == 0) boundary = .TRUE.
         
        IF (.NOT. boundary) THEN
         
           DO i = 1, 2
               
              ele_1 = edge_ele(i, iseg)
              
              is_1 = ele(ele_1)%verts(1)
              is_2 = ele(ele_1)%verts(2)
              is_3 = ele(ele_1)%verts(3)
              IF( SIZE(ele(ele_1)%NU_seg) == 4 ) THEN
                 is_4 = ele(ele_1)%verts(4)
              ENDIF
               
              cn_1 = cn(1, iseg) 
              cn_2 = cn(2, iseg)
               
              l_1 = ABS( (is_1 - cn_1)*(is_1 - cn_2) )
              l_2 = ABS( (is_2 - cn_1)*(is_2 - cn_2) )
              l_3 = ABS( (is_3 - cn_1)*(is_3 - cn_2) )
              IF( SIZE(ele(ele_1)%NU_seg) == 4 ) THEN
                 l_4 = ABS( (is_4 - cn_1)*(is_4 - cn_2) )
              ENDIF

              IF( SIZE(ele(ele_1)%NU_seg) == 3 ) THEN ! TRIANGLE
                 IF( (l_2 + l_3) == 0 ) ele(ele_1)%NU_seg(1) = iseg
                 IF( (l_1 + l_2) == 0 ) ele(ele_1)%NU_seg(3) = iseg
                 IF( (l_1 + l_3) == 0 ) ele(ele_1)%NU_seg(2) = iseg
              ELSE ! QUADRANGLE
                 IF( (l_3 + l_2) == 0 ) ele(ele_1)%NU_seg(1) = iseg
                 IF( (l_4 + l_3) == 0 ) ele(ele_1)%NU_seg(2) = iseg
                 IF( (l_1 + l_4) == 0 ) ele(ele_1)%NU_seg(3) = iseg
                 IF( (l_2 + l_1) == 0 ) ele(ele_1)%NU_seg(4) = iseg
              ENDIF                             
                          
            ENDDO
         
         ELSE
         
            ele_1 = edge_ele(1, iseg) + edge_ele(2, iseg)
            
            is_1 = ele(ele_1)%verts(1)
            is_2 = ele(ele_1)%verts(2)
            is_3 = ele(ele_1)%verts(3)
            IF( SIZE(ele(ele_1)%NU_seg) == 4 ) THEN
               is_4 = ele(ele_1)%verts(4)
            ENDIF
              
            cn_1 = cn(1, iseg)
            cn_2 = cn(2, iseg)
            
            l_1 = ABS( (is_1 - cn_1)*(is_1 - cn_2) )
            l_2 = ABS( (is_2 - cn_1)*(is_2 - cn_2) )
            l_3 = ABS( (is_3 - cn_1)*(is_3 - cn_2) )
            IF( SIZE(ele(ele_1)%NU_seg) == 4) THEN
               l_4 = ABS( (is_4 - cn_1)*(is_4 - cn_2) )   
            ENDIF            

            IF( SIZE(ele(ele_1)%NU_seg) == 3 ) THEN ! TRIANGLE
               IF( (l_2 + l_3) == 0 ) ele(ele_1)%NU_seg(1) = iseg
               IF( (l_1 + l_2) == 0 ) ele(ele_1)%NU_seg(3) = iseg
               IF( (l_1 + l_3) == 0 ) ele(ele_1)%NU_seg(2) = iseg
            ELSE ! QUADRANGLE
               IF( (l_3 + l_2) == 0 ) ele(ele_1)%NU_seg(1) = iseg
               IF( (l_4 + l_3) == 0 ) ele(ele_1)%NU_seg(2) = iseg
               IF( (l_1 + l_4) == 0 ) ele(ele_1)%NU_seg(3) = iseg
               IF( (l_2 + l_1) == 0 ) ele(ele_1)%NU_seg(4) = iseg
            ENDIF          

         ENDIF
             
      ENDDO

      ! n_ele(jt, i) : neighbour element to jt on the face i
      DO jt = 1, N_ele_mesh

         ALLOCATE( ele(jt)%n_ele(SIZE(ele(jt)%NU_seg)) )
                  
         DO i = 1, SIZE(ele(jt)%NU_seg)

            IF( edge_ele(1, ele(jt)%NU_seg(i)) == jt ) THEN
               ele(jt)%n_ele(i) = edge_ele(2, ele(jt)%NU_seg(i))
            ELSE
               ele(jt)%n_ele(i) = edge_ele(1, ele(jt)%NU_seg(i))
            ENDIF
            
         ENDDO 

      ENDDO
 
  END SUBROUTINE find_segments
  !===========================

  !========================================
  SUBROUTINE read_gmshMesh(unit, mesh_file)
  !========================================
  !
  ! Read 2D mesh file generated by Gmsh. 
  ! Node, triangle and quadriangles are read
  !

    IMPLICIT NONE

    INTEGER,           INTENT(IN) :: unit
    CHARACTER(len=64), INTENT(IN) :: mesh_file
    !==========================================

    INTEGER :: i, j, dummy, type_dummy
    INTEGER :: N_read, N_ele_tot, N_tri, N_qua
    INTEGER :: ierror
    !==========================================
              
    OPEN(UNIT, FILE = mesh_file, ACTION = 'READ', IOSTAT = ierror)
    IF (ierror /= 0) THEN
       WRITE(*,*) 'ERROR: Impossible to open mehs file ', TRIM(ADJUSTL(mesh_file))
       WRITE(*,*) 'STOP'     
       STOP
    ENDIF

    WRITE(*,*)
    WRITE(*,*) 'Getting mesh from ', TRIM(ADJUSTL(mesh_file)),' ...'
   
    ! Skip headers
    READ(unit, *); READ(unit, *); READ(unit, *)

    !----------------  
    ! Number of nodes
    !----------------------------------------------
    READ(unit, *)
    READ(unit, *) N_nodes
          
    ALLOCATE ( rr_nodes(N_dim, N_nodes) )
      
    DO j = 1, N_nodes     
       READ(unit, *) dummy, rr_nodes(:, j), dummy
    ENDDO
    
    READ(unit, *)
      
    WRITE(*,*) '   Found', N_nodes, 'nodes'   
    !----------------------------------------------
       
    ! Skip points
    READ(unit, *)
    READ(unit, *) N_ele_tot

    N_read = 0
    DO
       READ(unit, *) dummy, type_dummy, &
                     dummy, dummy, dummy, dummy
  
       IF (type_dummy /= GMSH_POINT) EXIT

       N_read = N_read + 1

    ENDDO 
    !
    BACKSPACE(unit)
    ! Skip edge    
    DO   
       READ(unit, *) dummy, type_dummy, & 
                     dummy, dummy, dummy, dummy, dummy
                       
       IF (type_dummy /= GMSH_LINE) EXIT

       N_read = N_read + 1

    ENDDO
    !
    BACKSPACE(unit)

    N_ele_mesh = N_ele_tot - N_read

    ALLOCATE( ele(N_ele_mesh) )

    !---------------
    ! Read triangles
    !----------------------------------------------------------
    N_tri = 0
    DO

       IF (N_read + N_tri == N_ele_tot) THEN
          READ(unit, *) 
          EXIT
       ENDIF

       READ(unit, *) dummy, type_dummy, dummy, dummy, dummy, &
                     dummy, dummy, dummy

       IF (type_dummy /= GMSH_TRIANGLE) EXIT
      
       N_tri = N_tri + 1

    ENDDO

    BACKSPACE(unit)

    IF( N_tri /= 0) THEN
        
        DO i = 1, N_tri
           BACKSPACE(unit)
        ENDDO

        DO i = 1, N_tri
           ALLOCATE( ele(i)%verts(3) )
    
           READ(unit, *) dummy, dummy, dummy, dummy, dummy, &
                          ele(i)%verts(3), ele(i)%verts(2), ele(i)%verts(1) ! square
                        ! ele(i)%verts(1), ele(i)%verts(2), ele(i)%verts(3) 
        ENDDO

    ENDIF
                
    WRITE(*,*) '   Found', N_tri, 'triangles'
    !----------------------------------------------------------

    !-----------------
    ! Read quadrangles
    !----------------------------------------------------------
    N_qua = N_ele_mesh - N_tri

    IF (N_qua /= 0) THEN

       DO i = N_tri+1, N_tri+N_qua
          ALLOCATE( ele(i)%verts(4) )

          READ(unit, *) dummy, dummy, dummy, dummy, dummy, &
                         ele(i)%verts(1), ele(i)%verts(4),  &
                       & ele(i)%verts(3), ele(i)%verts(2)  
       ENDDO
          
     ENDIF
                 
    WRITE(*,*) '   Found', N_qua, 'quadrangles'
    !----------------------------------------------------------

    CLOSE(unit) 

  END SUBROUTINE read_gmshMesh
  !===========================

  !=======================================
  SUBROUTINE read_MeshFBx(unit, mesh_file)
  !=======================================

    IMPLICIT NONE

    INTEGER,           INTENT(IN) :: unit
    CHARACTER(len=64), INTENT(IN) :: mesh_file
    !==========================================
    CHARACTER(len=64) :: Key
    INTEGER :: i, j, dummy
    INTEGER :: ierror
    !==========================================

    OPEN(UNIT, FILE = mesh_file, ACTION = 'READ', IOSTAT = ierror)
    IF (ierror /= 0) THEN
       WRITE(*,*) 'ERROR: Impossible to open mehs file ', TRIM(ADJUSTL(mesh_file))
       WRITE(*,*) 'STOP'     
       STOP
    ENDIF

    WRITE(*,*)
    WRITE(*,*) 'Getting mesh from ', TRIM(ADJUSTL(mesh_file)),' ...'

    READ(unit, *)
    READ(unit, *) Key

    DO

       SELECT CASE(TRIM(Key))

       CASE("Dimension")

          READ(unit, *) dummy !  N_dim
          READ(unit, *) Key

       CASE("Vertices")

          READ(unit, *) N_nodes

          ALLOCATE( rr_nodes(N_dim, N_nodes) )

          DO j = 1, N_nodes     
             READ(unit, *) rr_nodes(:, j), dummy
          ENDDO

          WRITE(*,*) '   Found', N_nodes, 'nodes'

          READ(unit, *) Key

       CASE("Edges")

          READ(unit, *) dummy ! skip edge

          READ(unit, *) Key

       CASE("Triangles")

             READ(unit, *) N_ele_mesh

             ALLOCATE( ele(N_ele_mesh) )

             DO i = 1, N_ele_mesh

                ALLOCATE( ele(i)%verts(3) )
                
                READ(unit, *) ele(i)%verts(:), dummy
 
             ENDDO

             WRITE(*,*) '   Found', N_ele_mesh, 'triangles'

             READ(unit, *) Key

       CASE("End")

          EXIT

       CASE DEFAULT

          WRITE(*,*) 'ERROR: unknown element'
          STOP

       END SELECT

    END DO

    CLOSE(unit)

  END SUBROUTINE read_MeshFBx
  !==========================  
  
END MODULE geometry
