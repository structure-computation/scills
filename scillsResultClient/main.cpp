#include <Soca/Com/SodaClient.h>
#include "ScillsResultUpdater.h"


int main( int argc, char **argv ) {
    // connection
    
    SodaClient sc( QHostAddress::Any, 8890 );
    if ( not sc.connected() ) return 1;

    //
    sc.reg_type( "ScillsAssemblyItem" );

    // attente
    while ( SodaClient::Event event = sc.event() ) {
        MP mp = event.mp();
        if ( mp.type() == "ScillsAssemblyItem" ) {
            ScillsResultUpdater scru;
            scru.exec( mp );
        } 
    }
}
