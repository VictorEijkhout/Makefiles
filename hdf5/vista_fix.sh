#!/usr/bin/env bash

if [ ! -f configure ] ; then
    echo "Please run this within the hdf5-1.14.4 directory."
    exit 1
fi

patch -p1 << '_EOF'
diff -ur a/m4/aclocal_fc.f90 b/m4/aclocal_fc.f90
--- a/m4/aclocal_fc.f90    2024-05-22 11:43:17.000000000 -0700
+++ b/m4/aclocal_fc.f90    2025-04-02 13:07:25.558777203 -0700
@@ -217,14 +217,21 @@
       ENDDO prec
 
       DO k = 1, num_rkinds
-         WRITE(stdout,'(I0)', ADVANCE='NO') real_kinds(k)
-         IF(k.NE.num_rkinds)THEN
-            WRITE(stdout,'(A)',ADVANCE='NO') ','
-         ELSE
-            WRITE(stdout,'()')
+         IF(real_kinds(k).GT.2)THEN
+            WRITE(stdout,'(I0)', ADVANCE='NO') real_kinds(k)
+            IF(k.NE.num_rkinds)THEN
+               WRITE(stdout,'(A)',ADVANCE='NO') ','
+            ELSE
+               WRITE(stdout,'()')
+            ENDIF
          ENDIF
       ENDDO
 
+      IF(real_kinds(1).EQ.2)THEN
+         num_rkinds = num_rkinds - 1
+      ENDIF
+
+
      WRITE(stdout,'(I0)') max_decimal_prec
      WRITE(stdout,'(I0)') num_ikinds
      WRITE(stdout,'(I0)') num_rkinds
_EOF
