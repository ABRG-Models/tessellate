/*
 * hexGeometry class
 * Author: John Brooke
 *
 * Date 2019/10
 *
 * Helper class for HexGrid and Dregion
 * to be included from Dregion
 * (maybe change to HexGrid later)
 */
#include<chrono>
using morph::HexGrid;
using namespace std;
using namespace morph;
using namespace std::chrono;
#define PI 3.1415926535897932

#define NE(hi) (this->Hgrid->d_ne[hi])
#define HAS_NE(hi) (this->Hgrid->d_ne[hi] == -1 ? false : true)

#define NW(hi) (this->Hgrid->d_nw[hi])
#define HAS_NW(hi) (this->Hgrid->d_nw[hi] == -1 ? false : true)

#define NNE(hi) (this->Hgrid->d_nne[hi])
#define HAS_NNE(hi) (this->Hgrid->d_nne[hi] == -1 ? false : true)

#define NNW(hi) (this->Hgrid->d_nnw[hi])
#define HAS_NNW(hi) (this->Hgrid->d_nnw[hi] == -1 ? false : true)

#define NSE(hi) (this->Hgrid->d_nse[hi])
#define HAS_NSE(hi) (this->Hgrid->d_nse[hi] == -1 ? false : true)

#define NSW(hi) (this->Hgrid->d_nsw[hi])
#define HAS_NSW(hi) (this->Hgrid->d_nsw[hi] == -1 ? false : true)

class hexGeometry
{
public:
    double tol = 0.0001;

    enum lineType  {VERTICAL, HORIZONTAL, SLANTED};

    const double INFINITE = 9999999.9999999;

    struct point {
        double first;
        double second;
    };


    struct  lineSegment {
        point start;
        point end;
    };

    struct line {
        double slope;
        double intercept;
    };

    /*
     * a directed line, angle is angle to x axis
     * intercept is cut on y axis if slope not infinite
     * intercept is cut on x axis if angle = +- pi/2
     */
    struct dLine {
        double angle;
        double intercept;
    };


    struct horzLine {
        double slope = 0;
        double intercept;
	};

    struct vertLine {
        double slope = 999.999;
        double xintercept;
    };

    struct cIntersect {
        bool intersect;
        point p;
    };

    struct segDist {
        bool inside = false;
        double dist = 999.999;
	};

    int value;


    //constructor
    hexGeometry() {
    }


    point pair2point(std::pair<double,double> inpair) {
        point result;
        result.first = inpair.first;
        result.second = inpair.second;
        return result;
    }

    std::pair<double, double> point2pair(point inPoint) {
        std::pair<double, double> result;
        result.first = inPoint.first;
        result.second = inPoint.second;
        return result;
    }

    point makePoint(double x, double y) {
        point result;
        result.first = x;
        result.second = y;
        return result;
    }


    /*!
     * Euclidean distance between two points
     */
    double distance(point a, point b) {
        double result;
        result = sqrt((a.first-b.first)*(a.first-b.first) + (a.second-b.second)*(a.second-b.second));
        return result;
    }

  //returns linetype from lineSegment
    lineType segmentLineType(lineSegment s) {
        lineType  result;
        double denominator = s.end.first - s.start.first;
        double numerator = s.end.first - s.start.first;
        if (denominator == 0)
        {
            return result = VERTICAL;
        }
        else if ( numerator == 0)
        {
            return result = HORIZONTAL;
        }
        else
        {
            return result = SLANTED;
        }

    }

    lineType lineLineType (line l) {
        lineType result;
        if (l.slope == INFINITE) {
            result = VERTICAL;
        }
        else if (l.slope == 0) {
            result = HORIZONTAL;
        }
        else {
            result = SLANTED;
        }
        return result;

    }
  // returns line from lineSegment, tests for verticality and horizontality
    line segment2line(lineSegment s)
	{
        line result;
        double denominator = s.end.first - s.start.first;
        double numerator = s.end.second - s.start.second;
        if (fabs(denominator) <= this->tol) {
           result.slope = INFINITE;
           result.intercept = s.end.first;
           return result;
        }
        else if (fabs(numerator) < this->tol) {
            result.slope = 0.0;
            result.intercept = s.end.second;
            return result;
        }
        else {
            result.slope = numerator / denominator;
            result.intercept =  (s.end.first * s.start.second - s.end.second * s.start.first)/denominator;
            return result;
        }
    }

    dLine segment2dLine (lineSegment s) {
        dLine result;
        double angle, intercept;
        double denominator = s.end.first - s.start.first;
        double numerator = s.end.second - s.start.second;
        if (numerator >= 0) {
            angle = atan2(numerator, denominator);
        }
        else {
            angle = atan2(numerator, denominator) + 2.0 * PI;
        }

        if (denominator == 0.0) {
            intercept = s.end.first;
        }
        else if (numerator >= 0.0) {
            intercept = s.end.second - s.end.first * tan(angle);
        }
        else {
            intercept = s.start.first + s.start.first * tan(angle);
        }
        result.angle = angle;
        result.intercept = intercept;
        return result;
    }


// returns a vertLine from a line segment
    vertLine segment2vert(lineSegment s)
    {
        vertLine result;
        result.xintercept = s.end.first;
        return result;
    }

 horzLine segment2horz(lineSegment s) {
   horzLine result;
   result.intercept = s.end.second;
   return result;
 }

	//returns lineSegment from two points
    lineSegment createLineSegment(point a, point b)
	{
        lineSegment result;
        result.start = a;
        result.end = b;
        return result;
	}

    //returns the point that is the start of a lineSegment
    point startSegment(lineSegment s) {
        return s.start;
    }

    //returns the point that is the end of a lineSegment
    point endSegment(lineSegment s) {
        return s.end;
    }

    //returns a signed angle between 2 lines
    double subtendLines( dLine start , dLine end) {
        return start.angle - end.angle;
    }

    //returns signed angle between two segments
    double subtendSegments(lineSegment s1, lineSegment s2) {
        dLine dL1 = segment2dLine(s1);
        dLine dL2 = segment2dLine(s2);
        return subtendLines(dL1,dL2);
    }

    //true if point is left of a lineSegment
    bool pointLeftSegment(point p, lineSegment s) {
        double phi, psi;
        lineSegment sp = createLineSegment(startSegment(s), p);
        psi = segment2dLine(s).angle;
        phi = segment2dLine(sp).angle;
        if (phi <= PI) {
            if ((phi <= psi) && (psi <= phi + PI)) {
                return true;
            }
            else {
                return false;
            }
        }
        else {
            if ((phi <= psi && psi <= 2*PI) || (psi <= phi - PI)) {
                return true;
            }
            else {
                return false;
            }
        }
    }


//returns a bool if point is on line or not
    bool pointOnLine(point p, line l) {
        bool result = false;
        double residual = p.second - p.first * l.slope - l.intercept;
        if (lineLineType(l) == VERTICAL && (fabs(p.first - l.intercept) < this->tol)) {
                result = true;
                return result;
        }
        if (residual < this->tol)
        {
           result = true;
        }
       return result;
    } //end of pointOnLine

    bool pointOnhorzLine(point p, horzLine h)
	{
        bool result;
        if ( (p.second - h.intercept) < this->tol)
        {
            result = 1;
        }
        else
        {
            result = 0;
        }
        return result;
	}

    bool pointOnvertLine(point p, vertLine h)
	{
        bool result;
        if ( fabs((p.first - h.xintercept)) < this->tol)
        {
            result = 1;
        }
        else
        {
            result = 0;
        }
        return result;
    }



//returns line through point p perpendicular to line l
    line perpPoint2Line (point p, line l) {
        line result;
        result.slope = - 1.0 / l.slope;
        result.intercept = p.second - result.slope * p.first;
        return result;
	}


// returns horizontal line through point perpendicular to vertLine
  horzLine perpPoint2vertLine( point p, vertLine h) {
    horzLine result;
	result.intercept = p.second;
	return result;
  }

// returns a vertLine through a point perpendicular to horzLine
  vertLine perpPoint2horzLine(point p, horzLine h)
  {
    vertLine result;
	result.xintercept = p.first;
	return result;
  }


//returns a structure giving the result of the intersection of two lines.
    cIntersect lineIntersect (line l1, line l2) {
        cIntersect result;
        lineType lt1 = lineLineType(l1);
        lineType lt2 = lineLineType(l2);
        if ((l1.slope == l2.slope) && (l1.intercept == l2.intercept)) {
            if (lt2 == VERTICAL && lt1 == VERTICAL) {
            //    std::cout << " lines vertical same intercept " << endl;
            }
            result.intersect = true;
            result.p.first = INFINITE;
            result.p.second = INFINITE;
            return result;
		}
        else if ((l1.slope == l2.slope) && (l1.intercept != l2.intercept)) {
            if (lt2 == VERTICAL && lt1 == VERTICAL) {
            //    std::cout << " lines vertical different intercept " << endl;
            }
            result.intersect = false;
            result.p.first = INFINITE;
            result.p.second = INFINITE;
            return result;
        }
        else if (lt1 == VERTICAL && lt2 != VERTICAL) {
            result.intersect = true;
            result.p.first = l1.intercept;
            // std::cout << " hex side is vertical " << l1.intercept << endl;
            result.p.second = l2.slope * l1.intercept + l2.intercept;
            return result;
        }
        else if (lt1 != VERTICAL && lt2 == VERTICAL) {
            result.intersect = true;
            result.p.first = l2.intercept;
            // std::cout << " hexagon side is vertical " << l2.intercept << endl;
            result.p.second = l1.slope * l2.intercept + l1.intercept;
            return result;
        }
        else {
            result.intersect = true;
            result.p.first = (l2.intercept - l1.intercept) / (l1.slope - l2.slope);
            result.p.second = (l1.slope*l2.intercept - l1.intercept*l2.slope) / (l1.slope - l2.slope);
            //std::cout << " neither side  is vertical " << l2.intercept << endl;
            return result;
        }
    }

    /*!
     * returns true if point is within the box defined by a line segment
     */
    bool pointInLineSegmentBox(point p, lineSegment s) {
        bool result = false;
        double xsign = (p.first - s.start.first) * (p.first - s.end.first);
        double ysign = (p.second - s.start.second) * (p.second - s.end.second);
        if (s.start.first > s.end.first - this->tol && s.start.first < s.end.first + this->tol) {
            xsign = 0.0;
            ysign = (p.second - s.start.second) * (p.second - s.end.second);
        }
        else if (s.start.second > s.end.second - this->tol && s.start.second < s.end.second + this->tol) {
            xsign = (p.first - s.start.first) * (p.first - s.end.first);
            ysign = 0.0;
        }
        else {
            xsign = (p.first - s.start.first) * (p.first - s.end.first);
            ysign = (p.second - s.start.second) * (p.second - s.end.second);
        }
        if ((xsign <= 0.0) && (ysign <= 0.0)) {
            result = true;
        }
        return result;
    }

    /*!
     * returns with bool member of cIntersect true if line segments intersect with their
     * point of intersection inside both bounding boxes. Point member given by the
     * line intersection.
     */
    cIntersect lineSegmentIntersect (lineSegment s1, lineSegment s2) {
        cIntersect result;
        line l1 =  segment2line(s1);
        line l2 = segment2line(s2);
        result = lineIntersect(l1, l2);
        //return result;
        if (result.intersect == true && pointOnLine(result.p,l1) && pointOnLine(result.p,l2)) {
            return result;
        }
        else {
            result.intersect = false;
            return result;
        }
    }


// returns a structure giving the result of the intersection of a line and a horzLine
  cIntersect horzIntersect (horzLine h, line l) {
    cIntersect result;
	result.intersect = 1;
	result.p.first = (h.intercept - l.intercept) / l.slope;
	result.p.second = h.intercept;
	return result;
 }

// returns a structure giving the result of the intersection of a line and a vertLine
  cIntersect vertIntersect (vertLine h, line l) {
    cIntersect result;
	result.intersect = 1;
	result.p.first = h.xintercept;
	result.p.second = h.xintercept * l.slope + l.intercept;
	return result;
  }

//gives the distance of a point from a line.
	 double pointLineDist (point p1, line l1) {
	   double result;
	   line l2 = perpPoint2Line(p1,l1);
	   cIntersect cI = lineIntersect(l1, l2);
	   point p2 = cI.p;
	   result = getdist(p1,p2);
	   return result;
	   }

//gives the distance of a point from a vertLine
  double pointVertDist(point p1, vertLine h1) {
    double result;
	result = fabs(p1.first - h1.xintercept);
	return result;
  }

 //gives the distance of a point from a horizontal line
   double pointHorzDist(point p1, horzLine h1) {
     double result;
	 result = fabs(p1.second - h1.intercept);
	 return result;
   }

     std::vector<lineSegment> hexSides(morph::Hex h) {
        vector<lineSegment> result;
        result.resize(6);
        double lr = static_cast<double> (h.getLR());
        double sr = static_cast<double> (h.getSR());
        point point1, point2;
        point1.first = h.x + sr; point1.second = h.y + lr/2.0;
        point2.first = h.x; point2.second = h.y + lr;
        result[0] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = h.x - sr; point2.second = h.y + lr/2.0;
        result[1] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = h.x - sr; point2.second = h.y - lr/2.0;
        result[2] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = h.x; point2.second = h.y - lr;
        result[3] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = h.x + sr; point2.second = h.y - lr/2.0;
        result[4] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = h.x + sr; point2.second = h.y + lr/2.0;
        result[5] = createLineSegment(point1, point2);
        //now return the vector of line segments
        return result;
    }


     std::vector<lineSegment> hexSides(point p, double d) {
        vector<lineSegment> result;
        result.resize(6);
        double lr = 2.0 * d / sqrt(3.0);
        double sr = d;
        point point1, point2;
        point1.first = p.first + sr; point1.second = p.second + lr/2.0;
        point2.first = p.first; point2.second = p.second + lr;
        result[0] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = p.first - sr; point2.second = p.second + lr/2.0;
        result[1] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = p.first - sr; point2.second = p.second - lr/2.0;
        result[2] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = p.first; point2.second = p.second - lr;
        result[3] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = p.first + sr; point2.second = p.second - lr/2.0;
        result[4] = createLineSegment(point1, point2);
        point1 = point2;
        point2.first = p.first + sr; point2.second = p.second + lr/2.0;
        result[5] = createLineSegment(point1, point2);
        //now return the vector of line segments
        return result;
    }

    bool pointInHexBox (point p, morph::Hex h) {
        bool result;
        double upperx = h.x + h.getSR() + this->tol;
        double lowerx = h.x - h.getSR() - this->tol;
        double uppery = h.y + h.getLR() + this->tol;
        double lowery = h.y - h.getLR() - this->tol;
        result = (p.first >= lowerx && p.first <= upperx && p.second >= lowery && p.second <= uppery);
        return result;
    }



    bool hexIntersectLineSegment(lineSegment s, point p, double d) {
        bool result = false;
        vector<lineSegment> hexSides = this->hexSides(p, d);
        for (int i=0; i<6; i++) {
            if (lineSegmentIntersect(hexSides[i], s).intersect) {
                result = true;
                break;
            }
        }
        return result;
    }

    bool hexIntersectLineSegment(lineSegment s, morph::Hex h) {
        bool result = false;
        point pt;
        pt.first = h.x;
        pt.second = h.y;
        vector<line> hexLines;
        for (int i=0; i<6; i++) {
            hexLines.push_back(this->segment2line(this->hexSides(h)[i]));
        }
        line l = this->segment2line(s);
        int intersected = 0;
        for (int i=0; i<6; i++) {
            pt = lineIntersect(hexLines[i], l).p;
            if (lineIntersect(hexLines[i],l).intersect && pointInHexBox(pt, h) && pointInLineSegmentBox(pt, s)) {
                intersected++;
                //cout << " line intersectd side i " << i <<  endl;
            }
        }
        if (intersected > 0) {
           //cout << "no of sides intersected " << intersected << endl;
           result = true;
        }
        return result;
    }
    vector <morph::BezCurvePath<float>> eqTriangleMesh(float d, vector<pair<double,double>>&  baryPoints, vector<morph::BezCurve<float>>& outer, pair<float, float> centroid = std::make_pair(0.0,0.0)) {
        cout << "just entering eqTriangleMesh" << endl;
        vector <morph::BezCurvePath<float>> result;
        result.resize(6);
        vector<pair<float,float>> p;
        pair<float,float> p0 = centroid;
        p.resize(6, std::make_pair(0.0f,0.0f));
        float pi6 = PI/6.0;
        for (int i=0; i<6; i++) {
           p[i] = std::make_pair(d * cos(pi6 + i * PI / 3.0), d * sin(pi6 + i * PI / 3.0));
           p[i].first = p[i].first + centroid.first;
           p[i].second =  p[i].second + centroid.second;
        }
        for (int i=0; i<6; i++) {
            morph::BezCurve<float> c0(p0, p[i]);
            morph::BezCurve<float> c1(p[i], p[(i+1)%6]);
            morph::BezCurve<float> c2(p[(i+1)%6], p0);
            result[i].addCurve(c0);
            result[i].addCurve(c1);
            result[i].addCurve(c2);
            double bp1 = (p0.first + p[i].first + p[(i+1)%6].first)/3.0;
            double bp2 = (p0.second + p[i].second + p[(i+1)%6].second)/3.0;
            baryPoints.push_back(std::make_pair(bp1,bp2));
            outer.push_back(c1);
        }
        return result;
    }

    vector<morph::BezCurvePath<float>> eqTriangleTess(double ds, vector<pair<double,double>>& centres, morph::BezCurvePath<float>& outerBound) {
        vector<morph::BezCurvePath<float>> result;
        vector<pair<double,double>> baryPoints;
        vector<morph::BezCurve<float>> outer;
        cout << "Just entering equTriangleTess" << endl;
        float pi3 = PI/3.0;
        float d = ds / sqrt(3.0);
        baryPoints.resize(0);
        result = eqTriangleMesh(d,baryPoints,outer);
        cout << "after eqTriangleMesh call 0 " << "baryPoints size " << baryPoints.size() << endl;
        for (int k=0; k<6; k++) {
            centres.push_back(baryPoints[k]);
        }
        for (int i=0; i<6; i++) {
            baryPoints.resize(0);
            outer.resize(0);
            pair<float,float> offset = std::make_pair(ds*cos(i*pi3), ds*sin(i*pi3));
            vector<morph::BezCurvePath<float>> basic = eqTriangleMesh(d, baryPoints, outer, offset);
            for (auto bp : basic) {
                result.push_back(bp);
            }
            for (int j=i; j<i+3; j++) {
                outerBound.addCurve(outer[(j+4)%6]);
            }
            for (int k=0; k<6; k++) {
                centres.push_back(baryPoints[k]);
            }
        }
        return result;
    }

    vector <morph::BezCurvePath<float>> isosTriangleMesh(double longSide, double shortSide, vector<vector<point>>&  vertices, point basePoint) {
        cout << "just entering isoTriangleMesh" << endl;
        vector <morph::BezCurvePath<float>> result;
        result.resize(2);
        vertices.resize(2);
        vector<pair<float,float>> p;
        p.resize(4, std::make_pair(0.0f,0.0f));
        float vertical = sqrt(longSide*longSide - shortSide*shortSide/4.0);
        p[0] = std::make_pair(-shortSide/2.0 + basePoint.first, basePoint.second);
        p[1] = std::make_pair(shortSide/2.0 + basePoint.first, basePoint.second);
        p[2] = std::make_pair(basePoint.first, vertical + basePoint.second);
        p[3] = std::make_pair(basePoint.first + shortSide, vertical + basePoint.second);
        morph::BezCurve<float> c0(p[0], p[1]);
        morph::BezCurve<float> c1(p[1], p[2]);
        morph::BezCurve<float> c2(p[2], p[0]);
        result[0].addCurve(c0);
        result[0].addCurve(c1);
        result[0].addCurve(c2);
        morph::BezCurve<float> c3(p[1], p[3]);
        morph::BezCurve<float> c4(p[3], p[2]);
        morph::BezCurve<float> c5(p[2], p[1]);
        result[1].addCurve(c3);
        result[1].addCurve(c4);
        result[1].addCurve(c5);
        vertices[0].push_back(pair2point(p[2]));
        vertices[0].push_back(pair2point(p[0]));
        vertices[0].push_back(pair2point(p[1]));
        vertices[1].push_back(pair2point(p[3]));
        vertices[1].push_back(pair2point(p[2]));
        vertices[1].push_back(pair2point(p[1]));
        cout << "vertices[0] size " << vertices[0].size() << endl;
        return result;
        //vertices[0][1] = pair2point(p[0]);
        //vertices[0][2] = pair2point(p[1]);
        //vertices[1][0] = pair2point(p[3]);
        //vertices[1][1] = pair2point(p[2]);
        //vertices[1][2] = pair2point(p[1]);
    }

    vector<morph::BezCurvePath<float>> isosTriangleTess(double ratio, const int rowX, const int rowY, vector<vector<point>>& vertices, morph::BezCurvePath<float>& outer) {
        vector<morph::BezCurvePath<float>> result;
        vector<vector<point>> baseVertices;
        vertices.resize(0);
        double spaceX = 0.16;
        double longSide = spaceX * ratio;
        double spaceY = sqrt(longSide*longSide - spaceX*spaceX/4.0);
        unsigned int seed;
        vector<pair<float,float>> p;
        p.resize(4);
        cout << "in isosTriangleTess rowX " << rowX << " rowY " << rowY << " spaceX " << spaceX << " spaceY " << spaceY << endl;
        double offset = 0;
        for (int j = -rowY; j < rowY + 1; j++) {
            offset = (1 - j) * (spaceX/2.0);
            for (int i = -rowX; i < rowX + 1; i++) {
                 baseVertices.resize(0);
                 point basePoint = this->makePoint(i*spaceX - offset, j*spaceY);
                 vector<morph::BezCurvePath<float>> basic = isosTriangleMesh(longSide, spaceX, baseVertices, basePoint);
                 cout <<  " in rowX rowY loop i " << i << " j " << j << " x " << i*spaceX << " j " << j*spaceY <<endl;
                 for (auto bp : basic) {
                     result.push_back(bp);
                 }
                 for (int k=0; k<2; k++) {
                     cout << "in k loop " << k << endl;
                     vertices.push_back(baseVertices[k]);
                     cout << " baseVertices size " << baseVertices.size();
                     cout << " baseVertices[k].size " << baseVertices[k].size() << endl;
                     cout << " k " << baseVertices[k][0].first <<  " , " << baseVertices[k][1].first << " , " << baseVertices[k][2].first << endl;
                 }
             }
        }
            p[0] = std::make_pair((-rowX - rowY/2 -1) *spaceX, -rowY*spaceY);
            cout << "p0.x " << p[0].first << " p0.y " << p[0].second << endl;
            p[1] = std::make_pair((rowX - rowY/2) *spaceX, -rowY*spaceY);
            cout << "p1.x " << p[1].first << " p1.y " << p[1].second << endl;
            p[2] = std::make_pair((rowX + rowY/2)*spaceX, rowY*spaceY);
            cout << "p2.x " << p[2].first << " p2.y " << p[2].second << endl;
            p[3] = std::make_pair((rowY/2 - rowX -1)*spaceX, rowY*spaceY);

        cout << "p3.x " << p[3].first << " p3.y " << p[3].second << endl;
        morph::BezCurve<float> c0(p[0], p[1]);
        morph::BezCurve<float> c1(p[1], p[2]);
        morph::BezCurve<float> c2(p[2], p[3]);
        morph::BezCurve<float> c3(p[3], p[0]);
        outer.addCurve(c0);
        outer.addCurve(c1);
        outer.addCurve(c2);
        outer.addCurve(c3);
        cout << "end of isosTriangleTess " << endl;
        return result;
    }

/*
 * returns a vector of points that represent the coordinates of the vertices of an isosceles array
 * If lPertub true then the vertices are randomly perturbed. This function works with triangle neigbours
 * which computes the topologica connectivity of the tessellation. If we set ratio = 1 we get a
 * tessellation of equal triangles but on an isosceles mesh. If we set pRatio non zero and lPerturb true then
 * we get scalene triangles meshed together
 */
    vector<point> isosVertices( double ratio, const int rowX, const int rowY, double pRatio, bool lPerturb = false) {
        vector<point> result;
        result.resize(0);
        double spaceX = 1.0 / (1.0 * (2 * rowX + 1));
        double longSide = spaceX * ratio;
        double spaceY = sqrt(longSide*longSide - spaceX*spaceX/4.0);
        unsigned int seed;
        chrono::milliseconds ms1 = chrono::duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch());
        seed = static_cast<unsigned int> (ms1.count());
        morph::RandUniform<double> ruf(seed);
        cout << "in isosPerturbTess rowX " << rowX << " rowY " << rowY << " spaceX " << spaceX << " spaceY " << spaceY << " pRatio " << pRatio << endl;
        int count=0;
        for (int j = -rowY; j < rowY + 2; j++) {
            double offset = j * (spaceX/2.0);
            for (int i = -rowX; i < rowX + 2; i++) {
                double radius = pRatio * spaceX * ruf.get() / 4.0;
                double angle = 2.0 * PI * ruf.get();
                double x = i*spaceX + offset;
                double y = j*spaceY;
                if (lPerturb) {
                    x = x + radius * cos(angle);
                    y = y + radius * sin(angle);
                }
                point basePoint = this->makePoint(x, y);
                result.push_back(basePoint);
                count++;
         cout<<" vertex " << count << " x " << x << " y " << y  << endl;
            }
        }
        return result;
    }

    /* This method takes the vector of vertex points and the topological arrangement of the triangular tessellation.
     * It returns a vector of BezCurvePaths representing the regions of the tessellation
     * and a BezCurvePath that gives the outer boundary which is a chain of line segments
     */

    vector<morph::BezCurvePath<float>> genTriangleTess(const int rowX, const int rowY, vector<point> vertices, vector<vector<int>> vIndices, morph::BezCurvePath<float>& outer) {
        vector<morph::BezCurvePath<float>> result;
        int numTriangles = (2*rowX+1) * (2*rowY+1) * 2;
        int stride = (2*rowX+1)*2;
        result.resize(numTriangles);
        vector<vector<int>> tessVertices;
        vector<morph::BezCurve<float>> edges;
        vector<std::pair<double,double>> triPoints;
        triPoints.resize(3);
        edges.resize(0);
        for (int j = -rowY; j < rowY + 1; j++) {
            int jrowY = j + rowY;
            for (int i = -rowX; i < rowY + 1; i++) {
                int idx = 2*(i+rowX);
                int index = idx + (jrowY)*stride;
                triPoints[0] = point2pair(vertices[vIndices[index][0]]);
                triPoints[1] = point2pair(vertices[vIndices[index][1]]);
                triPoints[2] = point2pair(vertices[vIndices[index][2]]);
         cout<<" region " << index << " v1 " << vIndices[index][0] << " v2 " << vIndices[index][1] << " v3 " << vIndices[index][2] << endl;
        cout << triPoints[0].first << " , " << triPoints[0].second << " | " << triPoints[1].first << " , " << triPoints[1].second << " | " << triPoints[2].first << " , " << triPoints[2].second << endl;

                morph::BezCurve<float> c0(triPoints[0],triPoints[1]);
                morph::BezCurve<float> c1(triPoints[1],triPoints[2]);
                morph::BezCurve<float> c2(triPoints[2],triPoints[0]);
                edges.push_back(c0);
                edges.push_back(c1);
                edges.push_back(c2);
                result[index].addCurve(c0);
                result[index].addCurve(c1);
                result[index].addCurve(c2);

                idx = 2*(i + rowX) + 1; //for down triangles
                index = idx + (jrowY)*stride;
                triPoints[0] = point2pair(vertices[vIndices[index][0]]);
                triPoints[1] = point2pair(vertices[vIndices[index][1]]);
                triPoints[2] = point2pair(vertices[vIndices[index][2]]);
         cout<<" region " << index << " v1 " << vIndices[index][0] << " v2 " << vIndices[index][1] << " v3 " << vIndices[index][2] << endl;
         cout << triPoints[0].first << " , " << triPoints[0].second << " | " << triPoints[1].first << " , " << triPoints[1].second << " | "<< triPoints[2].first << " , " << triPoints[2].second << endl;
                morph::BezCurve<float> c3(triPoints[0],triPoints[1]);
                morph::BezCurve<float> c4(triPoints[1],triPoints[2]);
                morph::BezCurve<float> c5(triPoints[2],triPoints[0]);
                edges.push_back(c3);
                edges.push_back(c4);
                edges.push_back(c5);
                result[index].addCurve(c3);
                result[index].addCurve(c4);
                result[index].addCurve(c5);
            }
        }

        for (int i = -rowX; i < rowX + 1; i++) {
            int idx = (i+rowX)*6 + 1;
            outer.addCurve(edges[idx]);
        }

        for (int j = -rowY; j < rowY +1; j++) {
            int idx = (j+rowY+1)*3*stride - 1;
            outer.addCurve(edges[idx]);
        }

        for (int i = rowX; i > -rowX - 1; i--) {
            int idx = (2*rowY)*3*stride + 6*(i+rowX+1) -3;
            outer.addCurve(edges[idx]);
        }

        for (int j = rowY; j > -rowY - 1; j--) {
            int idx = (j+rowY)*3*stride;
            outer.addCurve(edges[idx]);
        }
        return result;
    }




    vector<vector<int>> triangleNeighbors(const int rowX, const int rowY, vector<vector<int>>& vIndices) {
        vector<vector<int>> result;
        result.resize((2*rowX+1)*(2*rowY+1)*2);
        vIndices.resize((2*rowX+1)*(2*rowY+1)*2);
        cout << "rowX " <<  rowX << " rowY " << rowY << " result size " << result.size() << endl;
        int stride = 2*(2*rowX + 1);
        int vstride = 2*rowX + 2;
        int r1=0; int r2=0; int r3=0; int index=0;
        int v1=0; int v2=0; int v3=0;
        for (int j = -rowY; j < rowY + 1; j++) {
            int  jrowY = j + rowY;
            for (int i = -rowX; i < rowX + 1; i++) {
                int idx = 2*(i + rowX); //for up triangles
                r1 = idx - 1 +  (jrowY)*stride;
                r2 = idx + 1 + (j - 1 + rowY)*stride;
                r3 = idx + 1 + (jrowY)*stride;
                v1 = i + rowX + (jrowY + 1)*vstride;
                v2 = i + rowX + (jrowY)*vstride;
                v3 = i + rowX + 1 +  (jrowY)*vstride;
                if (j == -rowY) r2 = idx + (2*rowY)*stride + 1;
                if ((idx)%stride == 0) r1 = idx + stride -1;
                index = idx + (jrowY)*stride;
                cout << " i " << i << " j " << j << " region " << index << " r1 " << r1 << " r2 " << r2 << " r3 " << r3 << endl;
                cout << " i " << i << " j " << j << " region " << index << " v1 " << v1 << " v2 " << v2 << " v3 " << v3 << endl;
                result[index].push_back(r1);
                result[index].push_back(r2);
                result[index].push_back(r3);
                vIndices[index].push_back(v1);
                vIndices[index].push_back(v2);
                vIndices[index].push_back(v3);

                idx = 2*(i + rowX) + 1; //for down triangles
                r1 = idx - 1  + (j + 1 + rowY)*stride;
                r2 = idx - 1  + (jrowY)*stride;
                r3 = idx + 1 +  (jrowY)*stride;
                v1 = i + rowX + 1 + (jrowY+1)*vstride;
                v2 = i + rowX + (jrowY+1)*vstride;
                v3 = i + rowX + 1 + (jrowY)*vstride;
                if (j == rowY) r1 = idx -1;
                if (idx%stride == 2*(2*rowX+1) -1 ) r3 = (jrowY)*stride;
                index = idx + (jrowY)*stride;
                cout << "i " << i << " j " << j  << " region " << index << " r1 " << r1 << " r2 " << r2 << " r3 " << r3 << endl;
                cout << " i " << i << " j " << j << " region " << index << " v1 " << v1 << " v2 " << v2 << " v3 " << v3 << endl;
                result[index].push_back(r1);
                result[index].push_back(r2);
                result[index].push_back(r3);
                vIndices[index].push_back(v1);
                vIndices[index].push_back(v2);
                vIndices[index].push_back(v3);
            }
        }
        return result;
    }








//private:

//distance between two points
double  getdist(point p1, point p2) {
   double result;
   result = sqrt((p1.first - p2.first)*(p1.first - p2.first) + (p1.second-p2.second)*(p1.second-p2.second));
   return result;
   }


};


