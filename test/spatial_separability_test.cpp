#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <limits>
#define CMMCORE_DEBUG
#include <cmmcore/separability.h>

#include "cmmcore/utils.h"


using namespace cmmcore;
// Define a type for vectors
//std::vector<VectorND> set1 = {{-0.68271699086037385, -0.010340953914014683, 1.1541607007989305}, {-0.71065471113404888, -1.1512722359687020, -0.45179581338470343}, {-0.0014001349289074483, 0.71399631486879667, -1.0494752045479743}, {1.1815330110032800, 0.0045881114891104247, 0.58401240070027061}}
//;    std::vector<VectorND> set2 = {{-0.21810387855415669, 0.42309162549054058, 1.0315960439412366}, {0.24467567122539435, 1.5249868301442469, -1.0864228697150908}, {-0.43238232817874955, -2.8304117919050906, -0.9139107868649391}, {0.57714410932131432, -0.80022408019651348, 0.7917334052715711}}
//;    std::vector<VectorND> set3={{-2.7255170206814547, 1.0708885903017036, 1.0315960439412366}, {-2.2627374709019037, 2.1727837949554099, -1.0864228697150908}, {-2.9397954703060476, -2.1826148270939276, -0.9139107868649391}, {-1.9302690328059837, -0.15242711538535048, 0.7917334052715711}}
//;    // Run Task 1
//
// Task 1: Separating Planes
inline bool taskSeparatingPlanes(const std::vector<vec3>& set1, const std::vector<vec3>& set2) {

    const bool res=SAT3D(set1, set2);
    res?printf("Separable\n"):printf("Not Separable\n");
    return res;
}


int main() {
    // Example points for Task 1
  //  std::vector<lp::Vector> set1 = {
  //      {1.0, 2.0, 3.0},
  //      {1.5, 2.5, 3.5},
  //      {0.5, 1.5, 2.5}
  //  };
  //  std::vector<lp::Vector> set2 = {
  //      {-1.0, -2.0, -3.0},
  //      {-1.5, -2.5, -3.5},
  //      {-0.5, -1.5, -2.5}
  //  };
    std::vector<vec3> set1 = {
        {-0.68271699086037452, -0.010340953914014683, 1.1541607007989310},
        {-0.68970142092879327, -0.29557377442768645, 0.75267157225302239},
        {-0.69668585099721192, -0.58080659494135822, 0.35118244370711371},
        {-0.70367028106563079, -0.86603941545503016, -0.050306684838794857},
        {-0.71065471113404954, -1.1512722359687020, -0.45179581338470332}, {-0.51238777687750781, 0.17074336328168821, 0.60325172446220454}, {-0.44369277780806021, -0.087519264814796008, 0.40422785338078848}, {-0.37499777873861251, -0.34578189291128025, 0.20520398229937239}, {-0.3063027796691648, -0.60404452100776451, 0.0061801112179562939}, {-0.23760778059971732, -0.86230714910424877, -0.19284375986345975}, {-0.34205856289464132, 0.3518276804773911, 0.052342748125478344}, {-0.19768413468732729, 0.12053524479809441, 0.055784134508554677}, {-0.053309706480013243, -0.1107571908812023, 0.059225520891630981}, {0.09106472172730079, -0.34204962656049909, 0.062666907274707362}, {0.23543914993461479, -0.57334206223979578, 0.066108293657783701}, {-0.17172934891177474, 0.53291199767309394, -0.4985662282112478}, {0.048324508433405598, 0.32858975441098476, -0.29265958436367906}, {0.26837836577858598, 0.12426751114887558, -0.086752940516110336}, {0.48843222312376638, -0.08005473211323362, 0.1191537033314584}, {0.70848608046894679, -0.28437697537534273, 0.3250603471790271}, {-0.0014001349289081144, 0.71399631486879689, -1.0494752045479743}, {0.29433315155413864, 0.53664426402387533, -0.64110330323591302}, {0.59006643803718539, 0.3592922131789536, -0.23273140192385169}, {0.88579972452023237, 0.18194016233403199, 0.17564049938820958}, {1.1815330110032791, 0.0045881114891104247, 0.58401240070027072}};

    std::vector<vec3> set2 = {{-0.21810387855415669, 0.42309162549054058, 1.0315960439412366},
        {-0.10240899110926893, 0.69856542665396715, 0.50209131552715469},
        {0.01328589633561883, 0.9740392278173936, -0.027413412886927165},
        {0.12898078378050662, 1.2495130289808203, -0.55691814130100914},
        {0.24467567122539435, 1.5249868301442469, -1.0864228697150908},
        {-0.27167349096030491, -0.39028422885836728, 0.54521933623969254}, {-0.12180692303288509, -0.056792146004011328, 0.25469355193766313}, {0.028059644894534719, 0.27669993685034461, -0.035832232364366362}, {0.17792621282195453, 0.61019201970470061, -0.3263580166663958}, {0.32779278074937435, 0.94368410255905655, -0.61688380096842521}, {-0.32524310336645312, -1.2036600832072750, 0.058842628538148756}, {-0.14120485495650126, -0.81214971866198959, 0.0072957883481715611}, {0.042833393453450608, -0.42063935411670411, -0.044251051841805576}, {0.2268716418634025, -0.029128989571418629, -0.09579789203178278}, {0.41090989027335434, 0.36238137497386669, -0.14734473222175987}, {-0.37881271577260128, -2.0170359375561824, -0.42753407916339503}, {-0.16060278688011739, -1.5675072913199675, -0.24010197524131985}, {0.05760714201236649, -1.1179786450837528, -0.052669871319244742}, {0.27581707090485041, -0.66844999884753786, 0.13476223260283032}, {0.49402699979733428, -0.2189213526113232, 0.32219433652490548}, {-0.43238232817874955, -2.8304117919050906, -0.9139107868649391}, {-0.18000071880373361, -2.3228648639779461, -0.4874997388308116}, {0.072380890571282386, -1.8153179360508018, -0.061088690796683984}, {0.32476249994629841, -1.3077710081236575, 0.36532235723744366}, {0.57714410932131432, -0.80022408019651348, 0.7917334052715711}};
   std::vector<vec3> set3={{0.81411409601703122, 1.0708885903017036, 3.6645103245270798}, {0.92980898346191898, 1.3463623914651301, 3.1350055961129977}, {1.0455038709068067, 1.6218361926285567, 2.6055008676989160}, {1.1611987583516945, 1.8973099937919833, 2.0759961392848338}, {1.2768936457965823, 2.1727837949554099, 1.5464914108707521}, {0.760544483610883, 0.25751273595279572, 3.1781336168255354}, {0.91041105153830282, 0.59100481880715172, 2.8876078325235062}, {1.0602776194657226, 0.92449690166150766, 2.5970820482214765}, {1.2101441873931424, 1.2579889845158636, 2.3065562639194472}, {1.3600107553205623, 1.5914810673702195, 2.0160304796174175}, {0.70697487120473479, -0.55586311839611202, 2.6917569091239919}, {0.89101311961468666, -0.16435275385082659, 2.6402100689340147}, {1.0750513680246385, 0.22715761069445889, 2.5886632287440374}, {1.2590896164345904, 0.61866797523974437, 2.5371163885540602}, {1.4431278648445423, 1.0101783397850297, 2.4855695483640829}, {0.65340525879858657, -1.3692389727450194, 2.2053802014224479}, {0.87161518769107049, -0.91971032650880447, 2.3928123053445232}, {1.0898251165835544, -0.47018168027258977, 2.5802444092665984}, {1.3080350454760383, -0.020653034036374862, 2.7676765131886731}, {1.5262449743685222, 0.42887561219983983, 2.9551086171107483}, {0.59983564639243836, -2.1826148270939276, 1.719003493720904}, {0.85221725576745433, -1.6750678991667831, 2.1454145417550312}, {1.1045988651424703, -1.1675209712396388, 2.5718255897891589}, {1.3569804745174863, -0.65997404331249454, 2.9982366378232865}, {1.6093620838925022, -0.15242711538535048, 3.4246476858574142}}
;
   std::vector<vec3> set11 ={
       {0.90917781244586826, 1.4477811768454865, 4.3558873537850005},
       {0.90219338237744939, 1.1625483563318146, 3.9543982252390917},
       {0.89520895230903064, 0.8773155358181427, 3.5529090966931829},
       {0.88822452224061188, 0.59208271530447087, 3.1514199681472741},
       {0.88124009217219323, 0.30684989479079916, 2.7499308396013662},
       {1.0795070264287348, 1.6288654940411893, 3.8049783774482737}, {1.1482020254981824, 1.3706028659447047, 3.6059545063668579}, {1.2168970245676298, 1.1123402378482206, 3.4069306352854407}, {1.2855920236370775, 0.85407760975173641, 3.2079067642040249}, {1.3542870227065253, 0.59581498165525226, 3.0088828931226095}, {1.2498362404116015, 1.8099498112368924, 3.2540694011115479}, {1.3942106686189153, 1.5786573755575954, 3.2575107874946236}, {1.5385850968262293, 1.3473649398782988, 3.2609521738777003}, {1.6829595250335434, 1.1160725041990021, 3.2643935602607765}, {1.8273339532408577, 0.88478006851970536, 3.2678349466438532}, {1.4201654543944680, 1.9910341284325950, 2.7031604247748215}, {1.6402193117396482, 1.7867118851704857, 2.9090670686223903}, {1.8602731690848284, 1.5823896419083765, 3.1149737124699586}, {2.0803270264300089, 1.3780673986462673, 3.3208803563175273}, {2.3003808837751896, 1.1737451553841582, 3.5267870001650965}, {1.5904946683773347, 2.1721184456282980, 2.1522514484380952}, {1.8862279548603815, 1.9947663947833765, 2.5606233497501565}, {2.1819612413434282, 1.8174143439384545, 2.9689952510622173}, {2.4776945278264750, 1.6400622930935329, 3.3773671523742790}, {2.7734278143095219, 1.4627102422486116, 3.7857390536863402}}
;
    // Run Task 1
    bool correctResults[4]={false, true,true,false};
    bool res1=taskSeparatingPlanes(set1, set2);
    bool res2=taskSeparatingPlanes(set1, set3);
    bool res3=taskSeparatingPlanes(set11, set2);
    bool res4=taskSeparatingPlanes(set11, set3);

    for (auto& val: set11) {
        val*=1000.;
    }
    for (auto& val: set3) {
        val*=1000.;
    }
    for (auto& val: set2) {
        val*=1000.;
    }

    // Example vectors for Task 2
    bool res5=taskSeparatingPlanes(set11, set2);
    bool res6=taskSeparatingPlanes(set11, set3);
    assert(correctResults[0]==res1);
    assert(correctResults[1]==res2);
    assert(correctResults[2]==res3);
    assert(correctResults[3]==res4);
    assert(res3==res5);
    assert(res4==res6);


    return 0;
}