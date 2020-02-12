#pragma once 

// Define some util macros which help to declare a set of channel signal functions.
// A channel function has the following signature: double Func (double time)
// All of these free functions are assigned a channel Id integer which is used
// to query the channel signal data from a centralized controller.
#define CONCAT(a, ...) _CONCAT(a, __VA_ARGS__)
#define _CONCAT(a, ...) a ## __VA_ARGS__

#define COMPL(b) _CONCAT(COMPL_, b)
#define COMPL_0 1
#define COMPL_1 0

#define IIF(c) _CONCAT(IIF_, c)
#define IIF_0(t, ...) __VA_ARGS__
#define IIF_1(t, ...) t
#define A() 1

#define CHECK_N(x, n, ...) n
#define CHECK(...) CHECK_N(__VA_ARGS__, 0,)
#define PROBE(x) x, 1,

#define NOT(x) CHECK(_CONCAT(NOT_, x))
#define NOT_0 PROBE(~)

#define BOOL(x) COMPL(NOT(x))
#define IF(c) IIF(BOOL(c))

#define DECREMENT(x) _CONCAT(DEC_, x)
#define DEC_0 0
#define DEC_1 0
#define DEC_2 1
#define DEC_3 2
#define DEC_4 3
#define DEC_5 4
#define DEC_6 5
#define DEC_7 6
#define DEC_8 7
#define DEC_9 8
#define DEC_10 9
#define DEC_11 10
#define DEC_12 11
#define DEC_13 12
#define DEC_14 13
#define DEC_15 14
#define DEC_16 15
#define DEC_17 16
#define DEC_18 17
#define DEC_19 18
#define DEC_20 19
#define DEC_21 20
#define DEC_22 21
#define DEC_23 22
#define DEC_24 23
#define DEC_25 24
#define DEC_26 25
#define DEC_27 26
#define DEC_28 27
#define DEC_29 28
#define DEC_30 29
#define DEC_31 30
#define DEC_32 31
#define DEC_33 32
#define DEC_34 33
#define DEC_35 34
#define DEC_36 35
#define DEC_37 36
#define DEC_38 37
#define DEC_39 38
#define DEC_40 39
#define DEC_41 40
#define DEC_42 41
#define DEC_43 42
#define DEC_44 43
#define DEC_45 44
#define DEC_46 45
#define DEC_47 46

// Recursive expansion
#define EVAL(...)  EVAL1(EVAL1(EVAL1(__VA_ARGS__)))
#define EVAL1(...) EVAL2(EVAL2(EVAL2(__VA_ARGS__)))
#define EVAL2(...) EVAL3(EVAL3(EVAL3(__VA_ARGS__)))
#define EVAL3(...) EVAL4(EVAL4(EVAL4(__VA_ARGS__)))
#define EVAL4(...) EVAL5(EVAL5(EVAL5(__VA_ARGS__)))
#define EVAL5(...) __VA_ARGS__

#define EMPTY()
#define DEFER(id) id EMPTY()
#define OBSTRUCT(...) __VA_ARGS__ DEFER(EMPTY)()
#define EXPAND(...) __VA_ARGS__
#define NOOP(...)
#define WHEN(c) IF(c)(EXPAND, NOOP)

// Repeat a macro upto count
#define REPEAT(count, macro, ...) \
    WHEN(count) \
    ( \
        OBSTRUCT(REPEAT_INDIRECT) () \
        ( \
            DECREMENT(count), macro, __VA_ARGS__ \
        ) \
        OBSTRUCT(macro) \
        ( \
            DECREMENT(count), __VA_ARGS__ \
        ) \
    )
#define REPEAT_INDIRECT() REPEAT

// Macro to define channel functions:
#define FN_NAME_GEN(i, _) _DriveChannel##i, 
#define REGISTER_DRIVE_CHANNEL(channelId, pulseChannelProvider, _) \
    double _DriveChannel##channelId(double time) {\
        const int m_channelId = channelId;\
        return GetPulseValue(pulseChannelProvider, m_channelId, time);\
    }

#define DECLARE_CHANNEL_FUNC_ARRAY(N) channelFunctionType *g_channelFnArray[N]  = { EVAL(REPEAT(DECREMENT(N), FN_NAME_GEN, ~)) CONCAT(_DriveChannel, DECREMENT(N)) }

// ======================================================
// This is the only macro that we need to use
// Params:  
// - pulseChannelProvider: the pulse channel controller
// - N: number of channels that we want to declare
// This will create one function for each channel,
// each has a unique ID to retrieve waveform value from 
// the pulse controller.
// ======================================================
#define REGISTER_N_DRIVE_CHANNELS(pulseChannelProvider, N) \
    EVAL(REPEAT(N, REGISTER_DRIVE_CHANNEL, pulseChannelProvider, ~));\
    DECLARE_CHANNEL_FUNC_ARRAY(N);
