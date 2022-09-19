#include "ina.hpp"

ina::ina() {
        state1_ina      = 0.0;
        state2_ina      = 0.0202919327386704;
        state3_ina      = 0.936490564532038;
        gating          = (state1_ina * state1_ina * state1_ina)  * state2_ina * state3_ina;
}

double ina::Update_INa(double dt, double v){
        am      = (v == -47.13) * ( 3.2 )
                + (v != -47.13) * ( 0.32 * (v + 47.13) / (1.0 - exp( -0.1 * (v + 47.13) ) ) );
        bm      = 0.08 * exp( -v / 11.0 );

        ah      = (v >= -40.0) * ( 0.0 )
                + (v < -40.0) * ( 0.135 * exp( -( v + 80.0 ) / 6.8 ) );
        bh      = (v >= -40.0) * ( 1.0 / ( 0.13 * ( 1.0 + exp( -(v + 10.66) / 11.1 ) ) ) )
                + (v < -40.0) * ((3.56 * exp( 0.079 * v) + 3.1e5 * exp(0.35 * v)));

        aj      = (v >= -40.0) * (0.0) 
                +(v < -40.0) * (( ( -127140 * exp(0.2444*v) - 3.474e-5 * exp(-0.04391 * v)) * (v + 37.78)) /
                (1.0 + exp( 0.311 * (v + 79.23) ) ));
        bj      = (v >= -40.0) * ((0.3 * exp(-2.535e-7*v)) / (1.0 + exp( -0.1 * (v + 32.0) )))
                + (v < -40.0) * ((0.1212 * exp( -0.01052 * v )) / (1.0 + exp( -0.1378 * (v + 40.14) )));

        // Rush-Larsen (RL) mothod
        state1_ina = ( am/(am + bm) ) + (state1_ina - ( am/(am + bm) ) ) * exp(-(dt) * (am + bm) );
        state2_ina = ( ah/(ah + bh) ) + (state2_ina - ( ah/(ah + bh) ) ) * exp(-(dt) * (ah + bh) );
        state3_ina = ( aj/(aj + bj) ) + (state3_ina - ( aj/(aj + bj) ) ) * exp(-(dt) * (aj + bj) );

        gating          = ( state1_ina * state1_ina * state1_ina )  * state2_ina * state3_ina;

        return gating;
}