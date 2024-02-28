//  ntt64.c
//  Copyright (c) 2023 Plover Signature Team. See LICENSE.

//  === 64-bit Number Theoretic Transform

#ifndef POLYR_Q32

#include <stddef.h>
#include <stdbool.h>

#include "polyr.h"
#include "mont64.h"

//  === Roots of unity constants

/*
    (g is the smallest full generator for both q1 and q2)

    n   = 2048
    q
    r   = 2^64 % q
    g   = Mod(15, q)
    h   = g^(znorder(g)/(2*n))

    bitrev(n,i) = \
      b = binary(n + (i%n));\
      sum(2^(i-2) * b[i] for i in range(2, length(b)))

    Generating the Montgomery ("r") scaled constants plover_w_64[511]:

    w = vector(2047,i,lift(r*h^bitrev(n,i)))

    (   the function polyr_fntt() evaluates the input, intepreted
        as a polynomial with the 0-degree coefficient first,
        at the roots of unity h^(2*bitrev(n,i)+1) for i= 0, 1, 2,.. )

    e = vector(2047,i,lift(h^(2 *bitrev(n,i-1)+1)));
    evalv(v,x) = sum(i=0,2047,x^i*v[i+1])
    fntt(v) = vector(2048,i,lift(evalv(v, Mod(e[i],q))))

*/

// file generated with params/concrete-param-ring.py

static const int64_t plover_w_64[2047] = {
	   842590368738,   1832179252479,   1197335448900,
	  1649014791980,    592342242177,   1729585179059,    365642527020,
	   186401286637,   1326181061678,    373086876684,   1042967305735,
	   354932232724,    873508983567,   1403382312917,   1432078243024,
	  1942887938334,   1944687790500,    154861686718,   1150321521197,
	   237227504286,   1593203050871,    643965630923,     41826845932,
	  1340737482923,   1538967253090,    910031279782,   1779506649328,
	  1193664079614,    862469382130,    246485986736,   1408283543412,
	  1796835741702,   1687038543268,   1877759520164,   1386359146149,
	   391150629287,   1706845976935,   1937132548979,   1126743021947,
	  1936786383718,    794540046106,   1524506357099,    403047739172,
	  1632031449789,    534387054169,    889950580449,    865934872047,
	  1535139234505,    112869575928,    460655819873,    866073046577,
	  1147353598909,   1037104321353,    380918458833,   1991758399710,
	  1225889377876,   1432831266765,   1141359815745,   1272925400984,
	  1261099481224,   1355383276160,   1687872411424,    753306323498,
	  1499962746911,    423060640561,    295531983616,   1868042913766,
	  1356639715756,   1603699349986,   1083518620036,     29021159369,
	  1044683125261,   1446497513685,   1056829222538,   1302908224199,
	   284336949413,    686199753971,   1061478433164,    470715879602,
	  1224647487224,   1584756518443,    109277899104,    596690197307,
	  1027331700539,   1391784439263,    242097674368,   1146675597474,
	  1282910894754,   1900196555933,    780931830309,   1669698977751,
	   153706530768,    629255004180,   1776173025342,   1305923806132,
	   736235935312,   1490345403648,   1708095489678,   1973690154610,
	   526367754210,    802375965501,   1616326481027,    617695399697,
	   944451126422,   1036860099165,    614743206454,    757179622261,
	   248149214656,    680190314098,   1335048747490,   1179626691807,
	   800331501674,   1751752263233,    824762651138,    226749702345,
	  1477722929759,   1564210765397,   1665467189387,   1045128973624,
	  1634831926779,   1418149547798,    675835501253,   1714397208393,
	   349855153307,   1703685037316,    408317830050,   1961834219547,
	   216604495507,    857948610402,   1688192142502,   1007042464133,
	   541037521272,   1735124961715,    861996097727,    284730703825,
	  1489987312762,    610830574361,   1189950498635,   1178767295842,
	   785922318877,   1634229584311,    119094449328,    638748063673,
	  1444038483740,    653154513306,    765456993771,     67620083428,
	  1623586664315,     93838309744,   1968409454686,   1559263192386,
	   731290070235,   1357925688020,   1760949886913,    368028315902,
	  1691814784685,    876653620231,    203374145673,   1404356262072,
	    22284419098,    397633142827,    546287744823,    620015183278,
	   793068457507,   1004365547649,    973927613240,    244672288989,
	   696119174956,    274929361190,   1326954003464,    256671602928,
	  1676340437365,   1177096518662,   1211737623007,    287840113164,
	  1665808392108,      5065021664,   1443739272639,   1430045087837,
	    79695717844,    585918726587,   1688114900063,   1491803221405,
	   313716943971,    199815637884,   1392480001408,    949202859670,
	  1429042154178,    844676347473,   1135848391241,    324035964848,
	   408271469593,   1053930247861,     79662460012,    307632537257,
	  1666734005891,    656295938069,   1117356797231,   1989777152153,
	  1531526424352,   1602061223470,   1015654944254,   1849826243227,
	  1204169406435,   1894399298861,   1058940654372,   1347255996866,
	   238183009672,   1088674494907,    181970054112,   1367122347315,
	  1329948826744,    698917644381,   1118918352123,    205627711125,
	   607314308535,   1025184718672,    342652209294,    592037046787,
	  1749475421414,    411163994990,    165652874184,    727618923118,
	   520994820154,    733260521099,    325415462135,   1219760430117,
	  1151679774832,   1898466632521,    599473569721,    939681492736,
	   319889864905,    287213329571,   1554491878357,   1345315498376,
	  1287392882479,   2000824959628,   1653941293340,    140766215686,
	   738120975280,   1216635912939,    485189226403,   1976104705904,
	  1583377732217,   1258992436140,   1176271203635,   1311220258008,
	  1723726144406,   1413907130475,    712527656002,    635411645486,
	  1234518515005,    744485892455,   1764164345621,   1225535982240,
	  1585170176222,   1384306458575,   1781918637039,   1563977334258,
	   316130676386,    671088945115,   1727753319355,    406433210192,
	    96170721007,   1743447671711,    149842643913,    812910160455,
	  1800815737768,   1147621179278,   1245592421078,    707075782399,
	  1536174612820,     82508772591,   1447217609795,    239838936880,
	  1861697417644,   1423651026928,    295147852174,   1741572235479,
	   274013111599,   1640740343350,    370025592298,   1716641879361,
	   132905803303,   1662374882769,   1501297370558,    697613241373,
	  1719360952671,   1998391453060,   1200095999898,    561116248885,
	   654128414580,   1885819311728,    787290172345,    919572782557,
	  1639204540427,   1275110544866,   1931820566962,    114173298944,
	   930677127909,   1871285817910,    441594154989,    915244659606,
	  1219629814382,   1812212998783,   1624806337702,   1782668648657,
	   448243598600,      2314291510,    580003197932,   1754297577176,
	   297939028281,    316303105080,   1274804182460,    151358384479,
	   387651018246,   1581813290379,    418326771912,    369950157027,
	   423523633586,   1086294176245,   1910056710943,     90170484286,
	  1093626186946,   1439664046315,   1508449135864,   1035092688377,
	   838456266031,   1064945058559,    980451692003,   1418520746438,
	   942154453924,   1394501152526,   1021933053713,   1918844688708,
	   822993035136,    963246676386,   1641603556173,     76003703273,
	   638117627500,    782157008819,   1749353077988,    762970492958,
	  1045830898509,   1943628712625,   1033173724798,   1865809319241,
	   993970924610,   1666519578846,   1309315312698,   1032018922089,
	  1824040537450,    119539962944,    297428519715,   1898800281807,
	   124004695073,    739603168871,   1149641948142,   1006232524710,
	   391771334821,    586089017220,    833855074333,     30320304525,
	   247485076217,   1270818478777,   1768598392637,   1738258133207,
	  1002293203126,    780960518620,      6225766943,    809941080214,
	  1957197235215,   1741171504841,   1974050125078,   1686789011916,
	  1915855163904,   1502461392190,     62474859475,   1654453004173,
	  1694440500673,   1104363685076,   1100392905040,    850000267846,
	  1643145882756,    260189948242,    736569130302,   1428484676077,
	   588461026873,    928364423426,   1241141341340,   1750297159400,
	   990615309418,   1066607871196,   1538437101550,   1679618939052,
	   252143110573,    716977544754,   1316273744860,    261499413399,
	   545390382229,    203054077357,   1168012016151,    296558574893,
	  1088249807796,   1697970806282,    704161521454,   1474820082444,
	   604314136810,    792863074648,    540593286042,    272863231516,
	   893253588922,    766685457402,    293042956492,   1967059744508,
	   113411769759,   1834787668060,   1264330889272,   1916473335540,
	  1684296830632,   1829817380667,    473010212815,    935860765233,
	   839356001047,   1752160443638,   1901943970286,   1671576886219,
	   289395517299,    774693166057,   1995693827048,    546546063401,
	  1173559866786,     81876192981,   1698779381219,    814603795323,
	   968938680888,    326107424444,   1275940052734,   1138493922115,
	   649318231639,    364936429445,    204827149955,    295418964961,
	  1914979495437,   1061880799240,    415788424145,    256798970328,
	   941122460386,   1107646019903,   1674847215531,    150657124784,
	  1288977190836,   1781368968293,   1796782211649,     73999219870,
	   431810812233,   1136533036172,   1097850871020,    280715404233,
	  1637207840919,   1207714317250,    828214597208,    254723292890,
	   711608742864,   1100160531495,    499506152083,   1718525319479,
	  1849632967313,   1407805298239,    734769202379,   1630015276037,
	  1836667340197,    979047772657,    976311488075,   1308441550922,
	   576235836609,    882479228849,   1329331366181,    343408003429,
	    35154733272,   1632947422345,   1041737550482,    477508334765,
	   217676397658,    195352357538,   1276985548287,   1668620839432,
	  1110447530978,   1858178347324,    563899716904,    539572058648,
	   138213226095,   1978638555390,   1946257202559,   1738738997617,
	   849595287898,    869237761334,    935830089040,   1953486624887,
	   968366054689,   1997019399166,    497858104026,    721132248112,
	    34160825272,   1803886600133,    345959076666,   1333816530595,
	  1246901287149,    940611012478,    264322616601,    752128849868,
	  1052358636614,    790864749002,   1152589710139,   1009597183524,
	   118247508339,    409668539494,    650920812670,   1905247932865,
	  1671985048874,   1553761575561,   1453261896091,   1578749493699,
	  1009188440126,   1938517214549,   1927460000540,   1762895923503,
	  1135283665255,    499929741588,   1787986848704,   1550674520636,
	  1522505575107,   1071314299174,    719964241792,    469210674176,
	  1999514130956,    710494968732,   1699041030538,    665327559857,
	   874448602817,    474154383066,   1150141150400,   1336150018037,
	   744436352763,   1389343141431,    779636499182,    655124877078,
	    86049613059,   1891316237641,    203863521642,    612709241030,
	   812181118097,   1486028069488,   1164390712456,    159020441832,
	  1453279675597,   1592589921582,    971916651314,   1216135238999,
	   272389402825,   1639681739302,   2001393785782,      9773069908,
	  1925958928403,   1450442975529,    870854842228,   1032263563518,
	  1233847382227,   1919183765639,   1685784335617,   1988648226782,
	   788067068784,    921550796764,    343282841201,    270703127709,
	  1697604865405,    433687857516,   1906127344308,   1396038756044,
	  1087378669403,   1345424973649,   1792668449279,   1272398749735,
	  1193981626764,    829907930786,    611026654926,    536882448334,
	  1214306034111,   1360193104412,    311680249662,    521235378342,
	  1304130994475,   1889197583901,   1656616469200,    366032318210,
	   983833512634,    156384092685,    743045666735,   1888698487249,
	  1658325933333,    983401144984,   1327557887857,    197981271580,
	  1335636421098,    418803051455,   1022349461665,    999777293447,
	  1096264107514,   1224968591930,   1756346883711,    964970624188,
	  1637369045679,    166282720084,   1135168558422,    180647764408,
	   262275674455,    502259529005,    889428869171,   1824967150556,
	   438917886210,    404540744211,   1707411774229,   1115085275969,
	   611851908682,   1096811443584,   1182554898962,   1779234765861,
	   659409571355,    844027916953,    395104743110,   1581753763178,
	   904130838223,    495122975020,   1661899679984,    120838818151,
	   932589046546,    616288873239,    392787764762,   1682076437154,
	  1283665315406,    253806541657,   1974091183139,   1930025142715,
	   664528754305,    864600362743,   1970888585918,    165919977293,
	  1651241533829,   1772440189473,    211287819109,   1789258726747,
	  1966301512272,   1574760278107,   1734574292638,   1549870709054,
	   792287119216,   1005635384747,     23049344988,    219320773811,
	   989036732339,   1101586031976,    686508745699,    748064515301,
	   870267069387,   1434983946310,   1481448982515,    271612597481,
	   954187859600,   1818744431028,    975124181989,   1376882955596,
	  1575936290128,   1253680598025,    899575224380,    285580501809,
	  1008023521297,    781149874018,   1232298299227,   1930832043690,
	   507471782910,    594669404060,   1785278839517,   1689887482260,
	  1058120671332,    992700241041,   1156766883801,   1382276223704,
	   149822664105,   1784056112582,    485724119330,    631993566552,
	   444354313089,    109575361035,    644390562060,   1664046615454,
	   316648295343,    349520664568,   1172516441168,    313749235918,
	  1004776186445,    974702745219,   1001693728354,    461877745936,
	   225474463211,    917253058220,   1405745226016,   1301462223317,
	   358202111393,    458491085840,   1443336020614,   1424837185773,
	  1917336812676,    374244632019,   1621594409117,    208912892778,
	  1979746004968,   1609433489187,   1657868164467,    431995275636,
	   955562829580,   1714354150829,   1884310525539,    245599925408,
	   340425916127,   1355449666284,   1910479062628,   1261323287565,
	   489321555130,   1453918971893,    512258873618,   1691325015039,
	   467612440745,   1395890452843,    970438269724,    284920454121,
	   710787952085,    866663457703,   1151338204399,    433981685128,
	   573201006871,    868973024241,   1931166763673,   1107051525161,
	  1952312116149,   1181597150195,    529629838807,   1945303718413,
	  1972379199666,   1224883297872,     38603793824,    928023720510,
	   808259942587,   1898342858220,   1932022221795,    843195739713,
	   160045900765,   1441269760115,    860211213374,     69797527881,
	   504415839562,    346273340589,   1144724986146,    658772116061,
	  1891513128528,   1924405192715,   1954338689619,    759105305886,
	   780233728364,   1730840844606,   1983657347980,   1649227677141,
	   115171189334,    777548347549,    714853978590,    684787895573,
	   185264003636,     41523667355,    759286235714,   1620712573329,
	  1768641402098,   1931678558583,   1896163164068,   1457428664730,
	  1142043532825,   1105736167599,    595711805063,   1064910000083,
	  1291652333635,    406507021036,    212693575973,   1982708950221,
	  1879204775603,    450298437836,    946648705235,    898610755252,
	   105489774660,   1282816021657,   1990013138140,   1676954451872,
	  1106575129238,   1734185860071,   1585623714867,   1011757767558,
	  1366888128224,   1698333162274,    949519302981,   1122734374796,
	   996244495856,   1955695169287,    177635034130,    295485491751,
	   169081297210,    138536837290,   1176730869474,     87422539096,
	  1050650227475,   1404232347816,   1503045857250,   1490431977716,
	  1693298545694,     77704972534,    340348342836,   1406315479999,
	   146471160411,   1081467750331,     36635920769,    218489834947,
	   231757196532,   1917227612224,     82315267157,   1983743431038,
	   263713265519,   1575660151049,   1039530208473,   1163248003593,
	  1788204236824,   1335275093417,    841382957625,   1978541384609,
	   423466533373,    386524199330,    224042501959,    172304416581,
	   202319276079,   1740535407810,    753774348640,    984981316311,
	   263310866638,    222779962108,   1527069156180,   1844344384458,
	  1416746363121,   1250203254920,   1840548984896,   1966394521149,
	   300326729803,    599494365630,   1637135504682,    381920181413,
	   890754499871,   1505347972838,    187476710244,    182103420184,
	  1949727633666,     23810231488,    238286009443,    835261025002,
	  1572692696596,   1289352778265,   1883462792657,    395406430139,
	   584274027027,   1931749366289,    277175249704,    588237511502,
	  1013550878888,   1085595250890,   1267070393152,   1663153060464,
	   439457988357,   1756540804432,    717383744656,   1380294864437,
	  1096894081957,    573126393627,   1206123642953,    628787302745,
	   710493353821,    905617593907,   1235851514795,   1997275211584,
	   565615234880,   1906412163628,    356132934244,   1591512113366,
	  1811499345325,    279881961832,   1105237983125,    721400028123,
	   956311554194,    727916533055,    943620380122,   1414236407601,
	  1623527562163,   1947822229024,   1551594195261,    220009846773,
	  1999359854948,   1148480769135,   1821118232636,    461386011628,
	  1236277170485,    333236554567,    365781395088,   1418241151126,
	  1327065701913,    265050896766,   1016097266990,   1010609541432,
	   512216151630,   1315168477507,   1287034236215,    739276871203,
	  1187516497069,    428346533533,   1331867193564,   1478891212226,
	   889450061044,    143628099845,     36583786104,    319080194743,
	   817041323367,    107604028235,   1565259488781,    595042640526,
	  1020060154866,    366588695854,    451871939935,   1810550861738,
	  1601608282986,   1939686923116,   1528366728914,     99412072509,
	  1292608958932,   1346115309043,   1676960277461,   1792950298523,
	   631058878217,    524637910843,   1691337175186,    828161307934,
	   191747544771,   1149998811906,    494694370779,    644233838447,
	   878928488933,   1386145961048,   1436403227345,     61291364632,
	  1455154759353,     20670291858,   1272451860308,    729139679489,
	   857303573087,    397989070583,    204861464396,     35341448178,
	   557962205825,   1382498201945,   1100911310959,   1712136208971,
	   683275133575,   1487011018810,   1904229382821,   1069336954656,
	   698532963156,    163604339589,    784302370093,   1633341423183,
	   502103820024,   1672407567030,   1797718766720,   1601962408231,
	   667373334040,   1268870315358,   1762733117870,    466935729228,
	  1037660751402,   1740303496345,   1339759823354,    156368981376,
	  1377724495252,    957718313517,   1002491872170,   1557267495762,
	  1086782288569,   1118678008111,   1716529936656,   1664578725648,
	   212275318211,   1684731400391,   1580310195161,    107808233227,
	    11915733908,   1687752220045,    303809326593,    835724286907,
	  1650678683218,    294560553823,   1020560093820,   1657385321554,
	  1438443640096,    373270230793,   1697902717700,   1380783003965,
	   483029984034,   1026724649120,    183903809596,   1570481059516,
	    45982095849,   1056675840552,   1164331497016,   1790906771828,
	  1780080271591,   1833490623815,    757960363380,    262174009777,
	  1680821602331,   1679613127184,   1335652439991,   1881660954737,
	   532349319309,   1817592767389,    336472334530,    920910081909,
	   788409016488,    342311811973,    118533544052,   1361047321159,
	  1927296052362,   1309626500389,     26667200795,    924957554967,
	  1990925565043,    862374482545,    488947924448,    472276603900,
	   814261072734,   1633912433206,   1701210847534,   1582299589390,
	   509444557743,   1932974856576,    890503505323,    522520978822,
	  1647821152721,    922150897950,    599771812985,     34811982895,
	  1280648067139,    894041668014,   1535433948102,   1214642101056,
	  1566221659719,   1421430451606,   1194891895276,     67501279746,
	  1148587624061,   1048683677066,   1238069891611,    431465300801,
	  1622146793918,    319084760453,   1103961740195,   1155690292871,
	    51051671316,     11855869746,    904079268437,   1724444138394,
	   884986058249,     44545322181,    678157790854,    446837214256,
	   234797880605,    410860048991,    902037451638,   1787826368595,
	   365741881093,    281752632688,   1817404858371,    810665975307,
	  1313508992605,     15038686613,   1156840711186,    663278625213,
	  1603767562612,   1479026816071,    901753798901,    563460790477,
	   909248764432,   1867422351263,   1843543458753,    277187556716,
	  1281681556347,   1497848843060,   1390551072111,   1140855448730,
	  1401640650904,    331379609583,    620160639481,   1804159609138,
	  2001871544156,   1717424309594,   1322206820544,    349241330958,
	  1695117045211,    426980271292,   1418694413782,   1264594413213,
	   336083727060,    642359595871,   1448708519584,   1267618244651,
	   594800186598,   1032451283672,   1690507149658,   1542248801968,
	   588653001259,   1981311779411,   1089229545903,   1485867820307,
	   458464653886,      9829747355,   1661504219326,    552684186911,
	   389782910418,   1842383881405,    917943800405,   1871393132977,
	   942105390737,   1001120348791,    149498028247,   1501635042230,
	    61541777617,   1573804033967,   1033904828907,    639431730233,
	   314752164477,    648716372618,    434085284714,    107983150395,
	  1642133990697,   1179590715765,   1267385415308,   1714190970897,
	  1945468380769,   1362123093163,    721299107550,    653279305151,
	   463616098686,    956496411151,    707789633896,   1310896206723,
	  1245163126066,   1173899927934,   1777976763055,    475056947623,
	   342585064737,    929952701282,     11278345232,   1888009748149,
	   802842638683,   1947783619522,    868325688300,    928597293074,
	  1625527094183,   1015109828904,   1138286710132,   1515926576939,
	    67384397498,    571388905756,   1084040302161,   1124961161889,
	   672025428793,   1318093531150,   1252795489229,    215158976661,
	    73687911585,     34493814022,   1876822087824,    817709227901,
	   778776196040,    936678766989,   1506608397724,   1770727312063,
	  1559130386086,    114577510179,   1862227695076,    360582022367,
	   542456186877,   1981053000779,   1206940893765,   1932485854155,
	   183515892350,   1140708690874,   1348830217347,    428544443124,
	  1447920671181,   1027313992416,    557383948414,    826728705195,
	  1958898122083,   1057929746010,   1803239184494,   1989412520614,
	   434908117779,    868021480068,   1684069670017,   1680522997808,
	  1251382709874,    680529475222,   1415798521490,   1510446988440,
	   850833047164,    954821732309,   1876583086031,     76683530890,
	   544900750155,    171369412425,   1012644749444,    941587447494,
	   800611663793,    630453252240,    198759952220,   1105080607219,
	    69933348228,    105250673580,    657558117641,    439388300704,
	  1991779359190,    240157921196,    454039656875,    142273663777,
	  1299610522017,    628494677345,   1957136836546,   1207592110964,
	  1848286913944,    699872869223,   1958400012008,    897493888073,
	  1197317915124,   1697576251150,    612217354889,   1465089505422,
	   840578108829,   1456343080311,   1002276101143,   1765927117071,
	  1126357105586,   1487753322393,   1427120723540,   1122760335412,
	  1324516384617,   1604551353618,     31711522564,   1774487404485,
	  1278992531435,    749893748016,    658710952466,   1078182280950,
	   541856359486,   1892085953867,   1475449316178,   1552776768933,
	  1051176021899,    146098786655,   1222561400942,     91682055311,
	   639305852344,   1516707303040,    444249060268,    826711275112,
	   313947406796,    442094997625,   1788580152412,    466467385840,
	  1307221170221,     71659797182,   1388847844139,    981100642342,
	  1964848856288,   1535869618586,     14139459879,    574019785967,
	  1542751010762,   1155418052575,   1467461982442,    875452529540,
	   112899220325,    864557838697,    870874140590,   1350193208453,
	  1637093750077,   1342248570572,   1568992484128,   1557466726029,
	     1266537057,     72494593903,   1839681727278,    966485962721,
	  1227939299794,    101895639268,     69365812712,    310268683721,
	   150881125620,   1846894594743,   1912317269746,   1696741115544,
	  1553612054873,    401405324193,   1603123502166,   1879771398967,
	  1327143289186,   1071640646570,    206386476087,    138302253884,
	  1323306102006,    141789652731,    626965858884,    715696121885,
	   270419441255,    432797451794,   1525309761965,    404265564071,
	   103326478860,   1885255185457,    604980899289,    934779160773,
	  1449082440328,    194870583488,   1251313655593,    368047759270,
	   989437336969,    929693002066,   1127654204646,   1356845518218,
	   766071686486,   1784698617914,    721442868761,   1291011001168,
	    38738761667,    901755643629,    823678768137,   1104236800182,
	   818052544946,    532749157617,    878429396879,    825126166960,
	  1081197254071,    956218090815,    184766223014,   1094086856075,
	   651150499486,    290077374175,    901057210846,    990580062543,
	  1067598585254,    780303587898,   1493290412202,   1468368743991,
	  1049242931882,   1111975590707,    666315650889,     92806034934,
	   674381520014,   1955716915271,     37011888814,    508691514984,
	    20372817623,   1264179920070,    875249723714,   1381332317050,
	  1384514525777,   1348253279577,   1254029916097,   1410627671312,
	   487482648455,   1434170423269,   1541116237866,   1324063035712,
	  1445255701097,   1601648859286,   1593324553850,    926141138159,
	   860653101352,   1212914859829,   1647336642455,   1286167981108,
	   378257713374,    593031169774,    797117958799,   1990807628614,
	  1169324299847,    812614280610,    790176598486,   1228372964868,
	   605249647529,   1628315335237,   1301524373229,   1110283287741,
	   349801826164,    297777716034,    735355619142,    396145973036,
	   385908681559,    168107506502,    944671267002,    684982057394,
	  1353438678765,    697711589540,    459232838514,    875558274102,
	  1835670284383,   1251971601584,   1707384918071,    917732402992,
	  1323244888997,    130005313689,   1708672893060,   1725052821559,
	  1394047884351,    237370845763,   1858746193612,    355978183675,
	  1907912400920,    178254177110,   1293463525462,    250710504177,
	  1300783415293,   1500968511708,   1823864985585,   1339152001705,
	  1671304301659,   1049498325798,    244146081251,    753662724761,
	  1523746393207,     76188764086,    451134031390,   1536499027002,
	  1616715730104,   1585422630101,   1498449676166,    441859038157,
	  1121486457853,   1573612111932,    658946091921,    193978252396,
	  1050015196160,    626691569303,   1601796016617,     33217361691,
	   312236055986,   1849424389079,    917508994115,    144247369617,
	  1388282260729,   1111358893811,    133919778837,   1622818032835,
	   808620427175,   1211035695954,     94861753592,   1927677294275,
	  1423095846682,     66627200393,   1033210290837,    335871260525,
	  1210479624785,   1425752626866,    741458438312,   1477244355953,
	  1814108188164,   1851296824070,   1615370932077,    304867373732,
	   866881600350,    929347710534,    293003073294,   1345669308721,
	   158180971590,   1899391711471,   1570136086170,   1018800152084,
	    44429299153,    676087425043,    277757920052,    972743785339,
	   966496135216,   1671796673267,    484223845493,   1301400311288,
	   863750156890,   1362166748399,    231843789949,   1156714563261,
	    29132289790,   1714760756597,   1133850759758,    184670668576,
	  1770144461610,    108412950395,   1969248058390,   1747342014125,
	   893094458496,   1849692155723,   1982850212893,   1149700072711,
	  1838548777904,    518331582296,    927709880394,   1836752047969,
	  1042690347341,    655410123404,     93985071462,   1944749926186,
	    14536033689,   1913132315421,     78580873217,   1626988324780,
	  1924583936285,   1144536968315,    287718824916,   1903152612312,
	  1607764235442,   1971885178349,    463447465451,    135911671558,
	  1087506775694,    831644531325,    847004743316,    275789864665,
	  1084326837336,   1311492379868,   1220110085221,     48847935530,
	   893608579033,   1292241284785,   1716854809983,   1435648816349,
	  1873984642078,    255739031692,    639984953496,    933218667481,
	  1834949901548,    322501921372,    959431506437,   1020833841587,
	  1467946108031,    511391353351,    490840368848,    721033270610,
	   983018485385,    211564602576,   1494703395627,   1981157892773,
	  1791474419043,     75802517175,    425558058135,    921798267772,
	  1606692824506,   1221236547133,    216437986249,     48372856246,
	   937944440846,   1519097900612,    893251331264,   1637387049873,
	   607928249534,     57735953722,    168721741490,   1004748897893,
	  1036148285643,    831970624041,   1810255526340,    578015074952,
	   278614190389,   1565366844204,    424277813556,    896661949293,
	   449442694053,   1064327988522,    190846606563,     84377384038,
	  1707342463656,    437312229064,   1923739652712,    718626787672,
	   134652093909,   1860263732570,    368993561483,   1409574963509,
	    27731911578,    865912477884,   1761404785606,     95808076413,
	  1529603593440,   1940456124892,   1293671417051,    891072058382,
	  1877298299482,    326655026742,   1055120155862,    574751085726,
	  1938610198072,    105478173772,    380646506469,   1521528857299,
	  1517262118977,   1556374250073,    362416896548,   1828681987463,
	   938322514095,    510340743413,    679464464477,    467963122176,
	   367361819500,   1657193230734,   1782486644414,   1242859306463,
	   590025941121,    369660821566,   1089257360928,   1776062127321,
	  1193971659143,     68237671291,    797105316665,   1306034126836,
	  1428828042306,   1491140159109,   1473322178780,    363226551285,
	   372615775466,   1990667553216,   1821885166222,    714341745220,
	    61445881532,   1664110179752,    461314902628,    416755662276,
	  1464120433933,   1992863580599,   1365834381185,    528993864750,
	  1739315457173,   1539020048784,   1396673035250,    354526181777,
	  1254739345668,   1341368288148,    745350096493,    132375886478,
	  1880450519824,    828019266484,    525004404662,   1041385458435,
	   887189458877,   1999831907953,    745489810334,   1389463966994,
	   201597501541,   1037018763502,    571886985173,   1761121132703,
	   386165183201,    397999645744,   1804562745822,    879384126199,
	   681879162362,   1977053299262,    856451678836,   1614674281942,
	  1490409878210,    177543017278,    365425974699,    530547431501,
	   656281724712,    617049959569,   1594887108680,    370246678576,
	   534062758966,   1798750434673,   1802610082399,   1955179983058,
	   185297820986,    323401707096,   1177634871140,    408661332277,
	  1280754160718,   1458714798663,    761219724877,    458819688249,
	   135819004627,   1369764438924,    546000996277,   1206383754057,
	   802638910534,    613126734084,    936640345797,    497995792713,
	  1220281737825,    850894800446,   1577642169497,    341850217199,
	  1178981503438,    749517064902,   1256801445076,    621690943079,
	  1799346793294,    934615032191,   1525281464091,     61980790284,
	   516251519876,   1597407019622,   1025589698331,   1410989975344,
	   925677103381,   1989290411179,    929600612189,   1838330569146,
	  1923708608677,   1366325519429,   1561648616695,    884308527490,
	  1961222454364,    831082244770,    447893122946,   1050663535425,
	   833801301740,    501847031239,    574610940668,    933387408959,
	  1846338011377,   1166044553075,    983570238043,   1319695129893,
	  1560294876568,   1293015153742,   1221977643050,   1948948822682,
	   753824514612,   1230043048993,   1726348039523,   1007206372942,
	   215913788460,    761972485536,     59142615373,    919039987155,
	   908400537441,    528845828224,    731331592555,   1202775870372,
	  1838938178866,   1914545267762,   1881313925999,    592863749748,
	   743905971013,    260974368124,   1444119988959,   1671064298536,
	   213255158206,    269058262582,    670335803907,   1706131247027,
	  1222037643198,    268691931510,     25616546073,   1588874150144,
	  1977113297477,   1922377012143,   1137950011333,   1421577338842,
	    71655298023,   1484515856880,    867526748313,    164497568630,
	   817199128783,   1791811373846,    794791905879,    317993612253,
	   813898690222,   1634363559448,   1282740622266,   1909046958873,
	  1450514376600,    544759012272,    277942929404,    600523841213,
	  1310928894339,   1438302362987,   1084155057589,    694529235525,
	   234916706509,   1645111883953,    667958851015,    799062754643,
	   517264711111,    186364393388,   1230787882080,    524877836059,
	  1054635775457,   1492379213307,    623315613617,   1752816527932,
	  1597419505246,    735330441011,   1814460112490,    878144045230,
	   320883666929,    828297968141,   1230570976099,    130365855917,
	   415346444724,    598391669522,    527155718737,   1221361686316,
	   606460038148,    563988605355,   1472396635096,    877989526360,
	   766155911144,   1922027943070,    137403989227,    814486905776,
	  1978091677788,   1496308548868,   1572169511781,    981927470071,
	   720475439287,   1339368311047,    710691274691,    933514787517,
	  1707843699350,   1822234399919,   1374097106054,   1041666247846,
	   175810501873,    759410311586,    432182938369,    871254790185,
	  1469437713821,   1794267344665,   1477499528116,   1979146168972,
	   784332447371,   1921116682682,     64020824661,    327104127414,
	  1661433830660,   1070125898763,   1934639316773,   1443636728048,
	   601390430292,   1160368004875,     18496000814,    761775561972,
	  1863157937651,   1711800582260,   1385412711125,      3524348978,
	  1449902228328,   1892682344560,   1216726841362,   1137468319211,
	  1724714975664,     32957682361,    467174214837,   1384459128205,
	   357647721771,     56805294802,   1492137240547,   1049807398007,
	  1818845626669,    981121718680,    546118696017,   1027200104461,
	   943614601675,    553827345647,    647867077334,    802626369568,
	   768127084175,    918246279868,    634544529318,   1154514024710,
	    80786447373,   1794305279461,   1327255497224,   1107247500149,
	  1452830240242,   1397234937553,   1002698088029,     82599595881,
	  1891230305658,    274151246448,    260622760108,   1852663157238,
	  1093110425052,    294498740731,    497833207455,    840560528646,
	   124713298720,    238601477518,    363260498124,   1623959424046,
	  1294620596600,   1694752343743,    211280008275,    123541778791,
	  1764548809688,     90716790942,   1812378136388,    441927653568,
	  1239371738309,    600389611216,    830178865658,   1126952693620,
	   672013160579,    891339533360,   1116601787947,   1861762517282,
	  1524295282130,    530717349369,    460385262644,    253055344918,
	   567688763713,    225202362536,    167476896066,    553835878130,
	   421927770492,   1697747873731,   1404037906447,   1838767500080,
	   149019611288,    880095867898,     34317805669,    625249631720,
	   208517667429,    180949506574,    781319066460,   1028124263282,
	  1029060227605,    392783441121,   1394508107680,   1519175519056,
	  1884779628626,    359971859499,    450104016057,    812089885220,
	  1625972451184,    630456569192,   1319298359246,    603037866308,
	  1258753889237,   1381533121914,    129385662869,   1049498136064,
	   532663762882,    519960552204,   1203437043174,   1834891931368,
	  1322142527126,    979460954123,    505126122482,    560319143822,
	   402200568166,   1685739373993,     10545124336,    123709335965,
	  1427156855978,    327616242736,    546799619324,   1397117663474,
	  1261537136390,     96010466266,   1819094356466,    924622893997,
	  1470888486719,    678836942225,    473684003473,   1075429907590,
	   128689862458,    518179702112,   1830631375139,    483756715367,
	  1064779343447,    536778519298,   1014713272484,   1872216011668,
	   739698216281,   1409863021909,    319193750327,    847700165212,
	  1373531147322,   1678919075354,   1522881484530,   1314674524437,
	  1530637323545,   1429131269560,    127576551422,    243956802628,
	   181319075675,    795805553543,   1443363471674,    773099892755,
	  1793422350185,    670283511418,    667786692419,    827904475657,
	   589426285614,   1095678515365,   1772204840261,    275662663282,
	   271645219823,    779007540155,    672934613560,   1263508137476,
	   989422663129,    633569162000,   1505709594907,   1539376821805,
	   863470840390,   1205758118547,    257605651078,   1647700468742,
	    38024513635,   1450070625391,   1611622495434,   1278926893429,
	  1714679373739,    362109704329,   1351998025448,   1351071593505,
	  1184896866655,   1648932740597,    108933481302,   1846860550696,
};

// end generated      


//  Forward NTT (negacyclic -- evaluate polynomial at factors of x^n+1).

void polyr_fntt(int64_t *v)
{
    size_t i, j, k;
    int64_t x, y, z;
    int64_t *p0, *p1, *p2;

    const int64_t *w = plover_w_64;

    for (k = 1, j = PLOVER_N >> 1; j > 0; k <<= 1, j >>= 1) {

        p0 = v;
        for (i = 0; i < k; i++) {
            z = *w++;
            p1 = p0 + j;
            p2 = p1 + j;

            while (p1 < p2) {
                x = *p0;
                y = *p1;
                y = mont64_mulq(y, z);
                *p0++ = mont64_add(x, y);
                *p1++ = mont64_sub(x, y);
            }
            p0 = p2;
        }
    }
}

//  Reverse NTT (negacyclic -- x^n+1), normalize by 1/(n*r).

void polyr_intt(int64_t *v)
{
    size_t i, j, k;
    int64_t x, y, z;
    int64_t *p0, *p1, *p2;

    const int64_t *w = &plover_w_64[PLOVER_N - 2];

    for (j = 1, k = PLOVER_N >> 1; k > 0; j <<= 1, k >>= 1) {

        p0 = v;

        for (i = 0; i < k; i++) {
            z = *w--;
            p1 = p0 + j;
            p2 = p1 + j;

            while (p1 < p2) {
                x = *p0;
                y = *p1;
                *p0++ = mont64_add(x, y);
                y = mont64_sub(y, x);
                *p1++ = mont64_mulq(y, z);
            }
            p0 = p2;
        }
    }

    //  normalization
    polyr_ntt_smul(v, v, MONT_NI);
}

//  Scalar multiplication, Montgomery reduction.

void polyr_ntt_smul(int64_t *r, const int64_t *a, int64_t c)
{
    size_t i;

    for (i = 0; i < PLOVER_N; i++) {
        r[i] = mont64_cadd(mont64_mulq(a[i], c), PLOVER_Q);
    }
}

//  Coefficient multiply:  r = a * b,  Montgomery reduction.

void polyr_ntt_cmul(int64_t *r, const int64_t *a, const int64_t *b)
{
    size_t i;

    for (i = 0; i < PLOVER_N; i++) {
        r[i] = mont64_cadd(mont64_mulq(a[i], b[i]), PLOVER_Q);
    }
}

//  Coefficient multiply and add:  r = a * b + c, Montgomery reduction.

void polyr_ntt_mula(int64_t *r, const int64_t *a, const int64_t *b,
                    const int64_t *c)
{
    size_t i;

    for (i = 0; i < PLOVER_N; i++) {
        r[i] = mont64_csub(mont64_cadd(mont64_mulq(a[i], b[i]), PLOVER_Q) + c[i],
                           PLOVER_Q);
    }
}

//  POLYR_Q32
#endif
