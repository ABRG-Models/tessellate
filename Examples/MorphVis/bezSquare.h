pair<float,float> v1 = make_pair (-0.75f, -0.75f);
pair<float,float> v2 = make_pair (0.75f, -0.75f);
pair<float,float> v3 = make_pair (0.75f, 0.75f);
pair<float,float> v4 = make_pair (-0.75f, 0.75f);
//pair<float,float> v5 = make_pair (-0.3f, 0.1f);
cout << "ater making pairs" << endl;;
morph::BezCurve<float> c1(v1,v2);
morph::BezCurve<float> c2(v2,v3);
morph::BezCurve<float> c3(v3,v4);
//morph::BezCurve<float> c4(v4,v5);
morph::BezCurve<float> c4(v4,v1);
cout << "after making BezCurves" << endl;
morph::BezCurvePath<float> bound;
cout << "after making BezCurvePath" << endl;
bound.addCurve(c1);
cout << "after adding curve" << endl;
bound.addCurve(c2);
cout << "after adding curve" << endl;
bound.addCurve(c3);
cout << "after adding curve" << endl;
bound.addCurve(c4);
cout << "after adding curve" << endl;
//bound.addCurve(c5);
//cout << "after adding curve" << endl;
