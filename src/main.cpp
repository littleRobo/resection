#include <iostream>
#include "photogrammetry.h"

using namespace std;

void test_resection(){
    IO io;
    io.f = 1150;
    io.x0 = 225-227.5;
    io.y0 = 228.5-225;

    vector<coorC> cocs;
    coorC coc;
    coc.x = 30.990000-227.5; coc.y = 228.5-399.510000; cocs.push_back(coc);
    coc.x = 223.000000-227.5; coc.y = 228.5-387.940000; cocs.push_back(coc);
    coc.x = 43.250000-227.5; coc.y = 228.5-403.630000; cocs.push_back(coc);
    coc.x = 228.000000-227.5; coc.y = 228.5-367.920000; cocs.push_back(coc);
    coc.x = 30.490000-227.5; coc.y = 228.5-55.590000; cocs.push_back(coc);
    coc.x = 39.990000-227.5; coc.y = 228.5-57.810000; cocs.push_back(coc);
    coc.x = 221.990000-227.5; coc.y = 228.5-68.010000; cocs.push_back(coc);
    coc.x = 237.000000-227.5; coc.y = 228.5-78.910000; cocs.push_back(coc);

    vector<coorW> cows;
    coorW cow;
    cow.X = 239742.790000; cow.Y = 1188861.500000; cow.Z = 66.580000; cows.push_back(cow);
    cow.X = 240254.930000; cow.Y = 1188894.570000; cow.Z = 64.630000; cows.push_back(cow);
    cow.X = 239775.190000; cow.Y = 1188851.910000; cow.Z = 66.460000; cows.push_back(cow);
    cow.X = 240269.530000; cow.Y = 1188948.730000; cow.Z = 65.500000; cows.push_back(cow);
    cow.X = 239745.750000; cow.Y = 1189769.780000; cow.Z = 82.330000; cows.push_back(cow);
    cow.X = 239771.280000; cow.Y = 1189764.020000; cow.Z = 82.560000; cows.push_back(cow);
    cow.X = 240249.410000; cow.Y = 1189740.850000; cow.Z = 78.630000; cows.push_back(cow);
    cow.X = 240288.860000; cow.Y = 1189710.630000; cow.Z = 76.820000; cows.push_back(cow);

    EO eo;
    photogrammetry pgm;
    int i = pgm.resection(cocs, cows, io, eo);

    if(i){
        cout << "phi = " << eo.phi << endl;
        cout << "omg = " << eo.omg << endl;
        cout << "kaf = " << eo.kaf << endl;
        cout << "Xs = " << eo.Xs << endl;
        cout << "Ys = " << eo.Ys << endl;
        cout << "Zs = " << eo.Zs << endl;
    }

    //truth
    //phi= -0.0511735
    //omg= -0.0211784
    //kaf= 0.0034949
    //XS=  240413.6
    //YS= 1189391.5
    //ZS=    3088.3
}

void test_resection2(){
    IO io;
    io.f = 1150;
    io.x0 = 225-229.5;
    io.y0 = 229.5-225;

    vector<coorC> cocs;
    coorC coc;

    coc.x = 219.00-229.5; coc.y =229.5- 400.00; cocs.push_back(coc);
    coc.x = 409.75-229.5; coc.y =229.5- 387.75; cocs.push_back(coc);
    coc.x = 231.00-229.5; coc.y =229.5- 404.00; cocs.push_back(coc);
    coc.x = 414.75-229.5; coc.y =229.5- 368.00; cocs.push_back(coc);
    coc.x = 221.00-229.5; coc.y =229.5-  56.00; cocs.push_back(coc);
    coc.x = 231.00-229.5; coc.y =229.5-  58.25; cocs.push_back(coc);
    coc.x = 414.00-229.5; coc.y =229.5-  68.25; cocs.push_back(coc);
    coc.x = 428.50-229.5; coc.y =229.5-  79.25; cocs.push_back(coc);


    vector<coorW> cows;
    coorW cow;

    cow.X = 239742.790; cow.Y = 1188861.500; cow.Z = 66.580; cows.push_back(cow);
    cow.X = 240254.930; cow.Y = 1188894.570; cow.Z = 64.630; cows.push_back(cow);
    cow.X = 239775.190; cow.Y = 1188851.910; cow.Z = 66.460; cows.push_back(cow);
    cow.X = 240269.530; cow.Y = 1188948.730; cow.Z = 65.500; cows.push_back(cow);
    cow.X = 239745.750; cow.Y = 1189769.780; cow.Z = 82.330; cows.push_back(cow);
    cow.X = 239771.280; cow.Y = 1189764.020; cow.Z = 82.560; cows.push_back(cow);
    cow.X = 240249.410; cow.Y = 1189740.850; cow.Z = 78.630; cows.push_back(cow);
    cow.X = 240288.860; cow.Y = 1189710.630; cow.Z = 76.820; cows.push_back(cow);

    EO eo;
    photogrammetry pgm;
    int i = pgm.resection(cocs, cows, io, eo);

    if(i){
        cout << "phi = " << eo.phi << endl;
        cout << "omg = " << eo.omg << endl;
        cout << "kaf = " << eo.kaf << endl;
        cout << "Xs = " << eo.Xs << endl;
        cout << "Ys = " << eo.Ys << endl;
        cout << "Zs = " << eo.Zs << endl;
    }

    //truth
    //phi= 0.0459680
    //omg= -0.0789097
    //kaf= 0.0037945
    //XS=  239619.7
    //YS= 1189568.3
    //ZS=    3075.5
}

void test_registration(){
    vector<coorC> cocs1, cocs2;
    coorC coc;
    coc.x = 1; coc.y = 8; cocs1.push_back(coc);
    coc.x = 5; coc.y = 1; cocs1.push_back(coc);
    coc.x = 7; coc.y = 4; cocs1.push_back(coc);

    //truth
    //a1,a2,a3,b1,b2,b3
    //1,2,3,4,5,6

    coc.x = 20; coc.y = 50; cocs2.push_back(coc);
    coc.x = 10; coc.y = 31; cocs2.push_back(coc);
    coc.x = 18; coc.y = 54; cocs2.push_back(coc);

    photogrammetry pgm;
    regisEle reg;
    pgm.registration(cocs1,cocs2,reg);

    cout << "a1 = " << reg.a1 << endl
         << "a2 = " << reg.a2 << endl
         << "a3 = " << reg.a3 << endl
         << "b1 = " << reg.b1 << endl
         << "b2 = " << reg.b2 << endl
         << "b3 = " << reg.b3 << endl;

}

int main()
{
//    test_resection();
//    test_resection2();

    test_registration();

    cout << "Done!" << endl;
    return 0;
}
