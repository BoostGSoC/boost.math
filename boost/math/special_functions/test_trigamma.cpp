//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp> // Boost.Test
#include <boost/test/unit_test.hpp>
#include <boost/math/concepts/real_concept.hpp>
//#include <boost/math/special_functions/bernoulli.hpp>
#include "polygamma.hpp"
#include <boost/multiprecision/mpfr.hpp>

#define SC_(x) BOOST_STRINGIZE(x)

template <class T>
void test()
{
   static const char* data[] =
   {
       /* From 1 to 50 inclusive: */
      /* N[Table[polygamma(1,x),{x,1,50,1}],200] */
      SC_( 1.6449340668482264364724151666460251892189499012067984377355582293700074704032008738336289006197587053040043189623371906796287246870050077879351029463308662768317333093677626050952510068721400547968116),
      SC_( 0.64493406684822643647241516664602518921894990120679843773555822937000747040320087383362890061975870530400431896233719067962872468700500778793510294633086627683173330936776260509525100687214005479681156),
      SC_( 0.39493406684822643647241516664602518921894990120679843773555822937000747040320087383362890061975870530400431896233719067962872468700500778793510294633086627683173330936776260509525100687214005479681156),
      SC_( 0.28382295573711532536130405553491407810783879009568732662444711825889635929208976272251778950864759419289320785122607956851761357589389667682399183521975516572062219825665149398413989576102894368570045),
      SC_( 0.22132295573711532536130405553491407810783879009568732662444711825889635929208976272251778950864759419289320785122607956851761357589389667682399183521975516572062219825665149398413989576102894368570045),
      SC_( 0.18132295573711532536130405553491407810783879009568732662444711825889635929208976272251778950864759419289320785122607956851761357589389667682399183521975516572062219825665149398413989576102894368570045),
      SC_( 0.15354517795933754758352627775713630033006101231790954884666934048111858151431198494474001173086981641511543007344830179073983579811611889904621405744197738794284442047887371620636211798325116590792267),
      SC_( 0.13313701469403142513454668592040160645250999190974628354054689150152674477961810739371960356760451029266645048161156709686228477770795563374009160846238555120815054292785330804309681186080218631608594),
      SC_( 0.11751201469403142513454668592040160645250999190974628354054689150152674477961810739371960356760451029266645048161156709686228477770795563374009160846238555120815054292785330804309681186080218631608594),
      SC_( 0.10516633568168574612220100690805592744016431289740060452820121248918106576727242838137392455525883128032077146926588808451660576536227662139441259611670653886247153058217429569741779951512317397040692),
      SC_( 0.095166335681685746122201006908055927440164312897400604528201212489181065767272428381373924555258831280320771469265888084516605765362276621394412596116706538862471530582174295697417799515123173970406923),
      SC_( 0.086901872871768390750300180461774935704627122814755976429027658770172801304462511026002023728812550288585234279183243456417432211643268356931602678761334638036025249590438758507335154887024000416687915),
      SC_( 0.079957428427323946305855736017330491260182678370311531984583214325728356860018066581557579284368105844140789834738799011972987767198823912487158234316890193591580805145994314062890710442579555972243470),
      SC_( 0.074040268664010336838400114715555343331188595530074845594050669947030132007947060664397815970758638388519488059590870017890147530512433379942779536092038122585663645382680704595435089140804408043249387),
      SC_( 0.068938227847683806226155216756371669861800840428034029267520057702132172824273591276642713929942311857907243161631686344420759775410392563616248923847140163401990175994925602554618762610192163145290204),
      SC_( 0.064493783403239361781710772311927225417356395983589584823075613257687728379829146832198269485497867413462798717187241899976315330965948119171804479402695718957545731550481158110174318165747718700845759),
      SC_( 0.060587533403239361781710772311927225417356395983589584823075613257687728379829146832198269485497867413462798717187241899976315330965948119171804479402695718957545731550481158110174318165747718700845759),
      SC_( 0.057127325790782614376866481654487779050574389063174359909580803569106413500936413268184428655048040423843421554557484114509187303284287219517825240648370459441974797294425794788375010207270210050326728),
      SC_( 0.054040906037696194623780061901401359297487969310087940156494383816019993747849993515098008901961620670757001801471064361422767550197867466431405487561950706355555044208006041701955257120850456963906975),
      SC_( 0.051270822935203119831536294588381968715770517786542233785303248081947971587185173570499670951823116515632348061858875995771798021111994890254120169002393919651953936174765044472038359613925249207674288),
      SC_( 0.048770822935203119831536294588381968715770517786542233785303248081947971587185173570499670951823116515632348061858875995771798021111994890254120169002393919651953936174765044472038359613925249207674288),
      SC_( 0.046503249239057995114983006606522558284931515518968537640178531528659989727774742731497403378126971390915794773877016585340959018844421194108995452449105937792543505335762776898342214489208695919692429),
      SC_( 0.044437133536578656272007799994952310351047217998307380615385143098907923612072263392654428171515401142981910476356355428316165630414669127993292973110262962585931935087828892600821553332183902531262677),
      SC_( 0.042546774368336690298472828350339833980536821022882049802530700754862554992034456209289588852044701710089660948946147488807659014157580281112385600709506818918635148698414903942976562783979743741092544),
      SC_( 0.040810663257225579187361717239228722869425709911770938691419589643751443880923345098178477740933590598978549837835036377696547903046469170001274489598395707807524037587303792831865451672868632629981433),
      SC_( 0.039210663257225579187361717239228722869425709911770938691419589643751443880923345098178477740933590598978549837835036377696547903046469170001274489598395707807524037587303792831865451672868632629981433),
      SC_( 0.037731373316397176820497811913784935887177189201711767093786453549076887667905593618888536912531223735073224394048054129175837843874871536865179815042182690056044747646475390465001546347424845647732912),
      SC_( 0.036359631203914323596903847579079860441361002644784469425748044769927367777644962617516794800048370511479260059342978683359651286947573868826771035892662799795413746274733277982148322753460510942657467),
      SC_( 0.035084120999832690943842623089283942074014063869274265344115391708702877981726595270578019289844288878826198834853182764992304348172063664745138382831438309999495378927794502471944241120807449718167671),
      SC_( 0.033895060357739944213759388844337449802908237472128010885137983860902640169598176721232002642995299580371977669573753514100508866602503617182712699121569106670125581068103658238888355270629090621853759),
      SC_( 0.032783949246628833102648277733226338691797126361016899774026872749791529058487065610120891531884188469260866558462642402989397755491392506071601588010457995559014469956992547127777244159517979510742648),
      SC_( 0.031743366520302090126581680438741427141328864134169865434796903967273318860776347608039726079230702517127671969492819302052873301797323827611664022974037600137578465794661641820805339893128801571096446),
      SC_( 0.030766804020302090126581680438741427141328864134169865434796903967273318860776347608039726079230702517127671969492819302052873301797323827611664022974037600137578465794661641820805339893128801571096446),
      SC_( 0.029848530374755717307481588611376872504046954124987128979333175776272400587130801235220625987403337962490390059483636565597409573606322909338018476601218500045751101240024359910796157156673337842905445),
      SC_( 0.028983478471641530456270515947017010912351452394883322750959473354127071867407617844217165779790881215085545768826197119230627566685907684424523666912637185166858367676010519080346330167053960680275687),
      SC_( 0.028167151941029285558311332273547623157249411578556792138714575394943398398019862742176349453260268970187586585152727731475525525869581153812278768953453511697470612573969702753815717922156001496602218),
      SC_( 0.027395547002757680620039727335276018218977806640285187200442970456671793459748257803904744514988664031915981646881122793203920587597976215540673830681848573425865674302364764482210779650551063224997279),
      SC_( 0.026665086812838031240930888766977990461490589693608781064577375131617008945504284100471581622366311950104440375880392333014000938218867376972375802924361356479189268166499169157155995136307089521564116),
      SC_( 0.025972566037214762542869946938723142816061226812722354471779591198099003405338079114321997134831685911323276940977345241601258555947399232928054473284472159803288991158188919849676770759575787582505945),
      SC_( 0.025315103841291028157597100127414793046172995386029389317275975156021422866219078456859801211097300638476465632627575353369831862982244729312012395703933040802631528962265185464403923948267437812617713),
      SC_( 0.024690103841291028157597100127414793046172995386029389317275975156021422866219078456859801211097300638476465632627575353369831862982244729312012395703933040802631528962265185464403923948267437812617713),
      SC_( 0.024095219843670564148078956165487368893882691995190602880631120902600839879901410402130473430014611762807221135304553342661919905813892557985421081010298299577170493864109325857027362377773684094592728),
      SC_( 0.023528326419634282968940634170022516286172941428297178844349941764278844415048802692379906536590575481628082813309088490054210155246999133949139901871976304112317886154358758963603326096594545772597263),
      SC_( 0.022987493536995018501661023569698016556589382747929412484155241926528709206828142876263086633940494356695686923638996548464161480287561600147084736917947099136655365873125659991185803111197033603857404),
      SC_( 0.022470964611375183790917221916805454573118308367764123227956894819090692677902523041552342832287601794712215849258831259207963133180123583618159117083236355335002473311142188916805637821940835256749966),
      SC_( 0.021977137450881356630423394756311627412624481207270296067463067658596865517408695881058515671793774634218388688765004098714135972686296423124331956589409194841175312817315028422978477328113674762922805),
      SC_( 0.021504547658820865137039651845158508319996881963413963364249457072585523362399244085217305841926099775995326306912452113837009318622024211404105113489220158924351116219961531258517229691062635065380272),
      SC_( 0.021051854132338293837809230840178879528688597671879332309473540368194396155518302482682222093623700500304968679026530882510617286028090304658971568898907800391078142023492540765081285825059918904221377),
      SC_( 0.020617826354560516060031453062401101750910819894101554531695762590416618377740524704904444315845922722527190901248753104732839508250312526881193791121130022613300364245714762987303508047282141126443599),
      SC_( 0.020201333226697125805970645065733046773817941926588018505040202407159642117848812918148925781901732801660885195292901376286358875180758174528007618692975087169735183071204142412542991595803590522528564)
   };

   static const unsigned index[] =
   {
      1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50
   };

   static const unsigned table_size = sizeof(index) / sizeof(index[0]);

   T tol = boost::math::tools::epsilon<T>() * 10;
   for(unsigned i = 1; i < table_size; ++i)
   {
      T x(index[i]);
      T b2n = boost::math::trigamma(x);
      BOOST_CHECK_CLOSE_FRACTION(b2n, static_cast<T>(data[i]), tol);
   }

    std::cout<<typeid(T).name()<<std::endl;
}


BOOST_AUTO_TEST_CASE( test_main )
{
   using namespace boost::multiprecision;
   test<number<mpfr_float_backend<50>, et_off> >();
   test<number<mpfr_float_backend<100>, et_off> >();
   test<number<mpfr_float_backend<200>, et_off> >();
}


