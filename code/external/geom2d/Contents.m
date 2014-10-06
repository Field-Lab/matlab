% Geometry 2D Toolbox
% Version 0.3 28-Aug-2007 .
%
% 	Library to handle and visualize geometric primitives such as points,
% 	lines, circles and ellipses, polygons...
%
%   The goal is to provide a low-level library for manipulating geometrical
%   primitives, making easier the development of more complex geometric
%   algorithms. 
%
%   Most functions works for planar shapes, but some ones have been
%   extended to 3D or to any dimension.
%
% Points:
% -------
%   points2d                  - description of functions operating on points
%   clipPoints                - clip a set of points by a box
%   centroid                  - compute centroid (center of mass) of a set of points
%   polarPoint                - create a point from polar coordinates (rho + theta)
%   angle2Points              - return horizontal angle between 2 points
%   angle3Points              - return oriented angle made by 3 points
%   angleSort                 - sort points of plane according to their angles
%   distancePoints            - compute distance between two points
%   minDistancePoints         - compute minimal distance between several points
%
% Vectors:
% --------
%   vectors2d                 - description of functions operating on plane vectors
%   vecnorm                   - compute norm of vector or of set of vectors
%   vectorAngle               - compute angle of a vector with horizontal axis
%   normalize                 - normalize a vector
%   isPerpendicular           - check orthogonality of two vectors
%   isParallel                - check parallelism of two vectors
%
% Line objects (lines, edges and rays):
% -------------------------------------
%   lines2d                   - description of functions operating on planar lines
%   edges2d                   - description of functions operating on planar edges
%   rays2d                    - description of functions operating on planar rays
%   createEdge                - create an edge between two points, or from a line
%   createLine                - create a line with various inputs.
%   createRay                 - create a ray (half-line)
%   medianLine                - create a median line between two points
%   cartesianLine             - create a line with cartesian coefficients
%   orthogonalLine            - create a line orthogonal to another one.
%   parallelLine              - create a line parallel to another one.
%   bisector                  - return the bisector of two lines, or 3 points
%   lineFit                   - least mean square line regression
%
% Operations on line objects:
% ---------------------------
%   clipEdge                  - clip an edge with a rectangular box
%   clipLine                  - clip a line with a box
%   invertLine                - return same line but with opposite orientation
%   edgeAngle                 - return angle of edge
%   edgeLength                - return length of an edge
%   lineAngle                 - return angle between lines
%   linePosition              - return position of a point on a line
%
% Points and lines:
% -----------------
%   intersectEdges            - return all intersections points of N edges in 2D
%   intersectLineEdge         - return intersection between a line and an edge
%   intersectLinePolygon      - get intersection pts between a line and a polygon
%   intersectLines            - return all intersection points of N lines in 2D
%   distancePointEdge         - compute distance between a point and an edge
%   distancePointLine         - compute distance between a point and a line
%   projPointOnLine           - return the projection of a point on a line
%   pointOnLine               - create a point on a line at a given distance from line origin
%   onEdge                    - test if a point belongs to an edge
%   onLine                    - test if a point belongs to a line
%   onRay                     - test if a point belongs to a ray
%   isLeftOriented            - check if a point is on the left side of a line
%
% Polygons:
% ---------
%   polygons2d                - description of functions operating on polygons
%   clipPolygon               - clip a polygon with a rectangular window
%   clipPolygonHP             - clip a polygon with a Half-plane defined by a directed line
%   polygonCentroid           - compute centroid (center of mass) of a polygon
%   polygonExpand             - 'expand' a polygon with a given distance
%   polygonArea               - compute the signed area of a polygon
%   polygonLength             - compute perimeter of a polygon
%   polygonNormalAngle        - compute normal angle at a vertex of the polygon
%   polygonBounds             - compute bounding box of a polygon
%   polygonContains           - test if a point is contained in a multiply connected polygon
%   splitPolygons             - convert a NaN separated polygon list to cell array of polygons
%   readPolygon               - read a polygon stored in a file
%   rectAsPolygon             - convert a (centered) rectangle into a series of points
%   steinerPoint              - compute steiner point (weighted centroid) of a polygon
%   steinerPolygon            - create a Steiner polygon from a set of vectors
%   supportFunction           - compute support function of a polygon
%   convexification           - compute convexification of a polygon
%   medialAxisConvex          - compute medial axis of a convex polygon
%
% Polylines:
% ----------
%   polylines2d               - description of functions operating on polylines
%   curveLength               - return length of a curve (a list of points)
%   curveCentroid             - compute centroid of a curve defined by a series of points
%   parametrize               - return a parametrization of a curve
%   curvature                 - estimate curvature of a curve defined by points
%   subCurve                  - extract a portion of a curve
%   cart2geod                 - convert cartesian coordinates to geodesic coord.
%   geod2cart                 - convert geodesic coordinates to cartesian coord.
%
% Circles and Ellipses:
% ---------------------
%   circles2d                 - description of functions operating on circles
%   createCircle              - create a circle from points
%   createDirectedCircle      - create a directed circle
%   circleAsPolygon           - convert a circle into a series of points
%   ellipseAsPolygon          - convert an ellipse into a series of points
%   circleArcAsCurve          - convert a circle arc into a series of points
%   enclosingCircle           - find the minimum circle enclosing a set of points.
%   onCircle                  - test if a point is located on a given circle.
%   inCircle                  - test if a point is located inside a given circle.
%
% Polynomial curves:
% ------------------
%   polynomialCurves2d        - description of functions operating on polynomial curves
%   polynomialCurveCentroid   - compute the centroid of a polynomial curve
%   polynomialCurveCurvature  - compute the local curvature of a polynomial curve
%   polynomialCurveCurvatures - compute curvatures of a polynomial revolution surface
%   polynomialCurveDerivative - compute derivative vector of a polynomial curve
%   polynomialCurveProjection - projection of a point on a polynomial curve
%   polynomialCurveFit        - fit a polynomial curve to a series of points
%   polynomialCurveSetFit     - fit a set of polynomial curves to a segmented image
%   polyfit2                  - polynomial approximation of a curve
%   polynomialCurveLength     - compute the length of a polynomial curve
%   polynomialCurveNormal     - compute the normal of a polynomial curve
%   polynomialCurvePoint      - compute point corresponding to a position
%   polynomialCurvePosition   - compute position on a curve for a given length
%   polynomialDerivate        - derivate a polynomial
%
% Other shapes:
% -------------
%   squareGrid                - generate equally spaces points in plane.
%   hexagonalGrid             - generate hexagonal grid of points in the plane.
%   triangleGrid              - generate triangular grid of points in the plane.
%   crackPattern              - create a (bounded) crack pattern tessellation
%   crackPattern2             - create a (bounded) crack pattern tessellation
%
% Geometric transforms:
% ---------------------
%   transforms2d              - description of functions operating on transforms
%   homothecy                 - create a homothecy as an affine transform
%   lineSymmetry              - create line symmetry as 2D affine transform
%   rotation                  - return 3*3 matrix of a rotation
%   translation               - return 3*3 matrix of a translation
%   scaling                   - return 3*3 matrix of a scale in 2 dimensions
%   transformPoint            - transform a point with an affine transform
%   transformVector           - transform a vector with an affine transform
%   transformEdge             - transform an edge with an affine transform
%   transformLine             - transform a line with an affine transform
%
% Utility objects:
% ----------------
%   formatAngle               - Ensure an angle value is comprised between 0 and 2*PI
%   boxes2d                   - conventions for using bounding boxes
%
% Drawing functions:
% ------------------
%   drawArrow                 - draw an arrow on the current axis
%   drawCenteredEdge          - draw an edge centered on a point
%   drawCircle                - draw a circle on the current axis
%   drawCircleArc             - draw a circle arc on the current axis
%   drawCurve                 - draw a curve specified by a list of points
%   drawEdge                  - draw the edge given by 2 points
%   drawEllipse               - draw an ellipse on the current axis
%   drawEllipseArc            - draw an ellipse on the current axis
%   drawParabola              - draw a parabola on the current axis
%   drawLabels                - draw labels at specified positions
%   drawLine                  - draw the line on the current axis
%   drawPoint                 - draw the point on the axis.
%   drawPolygon               - draw a polygon specified by a list of points
%   drawRay                   - draw a ray on the current axis
%   drawRect                  - draw rectangle on the current axis
%   drawRect2                 - draw centered rectangle on the current axis
%   drawShape                 - draw various types of shapes (circles, polygons ...)
%   fillPolygon               - fill a polygon specified by a list of points
%
%
%   Credits:
%   * function 'enclosingCircle' rewritten from a file from Yazan Ahed
%       (yash78@gmail.com), available on Matlab File Exchange
%
%   -----
%
%   author : David Legland
%   INRA URPOI (Nantes) & MIA (Jouy-en-Josas)
%   david.legland@nantes.inra.fr
%   created the 07/11/2005.
%   Licensed under the terms of the LGPL, see the file "license.txt'


help('Contents');

%   Requires further development:
%   inertiaEllipse            - inertia ellipse of a set of points

%   Deprecated functions:
%   createMedian              - create a median line
%   minDistance               - compute minimum distance between a point and a set of points
%   clipLineRect              - clip a line with a polygon
%   triMeshInterp             - interpolate a tria

