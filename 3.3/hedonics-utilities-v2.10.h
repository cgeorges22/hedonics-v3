#ifndef UTILITIES_H
#define UTILITIES_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <iostream>
#include <fstream>

//Function Declarations

//initialize firms 
void initializeFirms(firm firms[], const int firmNum, bool twoClasses, std::ofstream &testoutput, bool anotation, double probComplements, double probInnovate, double probIntegerConstrained, double intPriceDiscount, double wage, double salary, double markup); //int round removed 6-14-15

//random number, uniform distribution on [0,1]
double randNumUni();

//N(0,1) random number
double randNumNorm();

//exogenous product innovation
void updateHedonicQual(firm firms[], const int firmNum, std::ofstream &testoutput, int round, double probInnovate, double probUpTick, double probDownTick, bool anotation);

//endogenous product innovation -- update firm's hedonic qualities conditional on recent R&D investment 
void doEndogProductInnovation(firm firms[], const int firmNum, std::ofstream &testoutput, int round, double probInnovate, double probUpdateZeroHQElement, double probUpTick, double probDownTick, bool multiplicativeMutation, double incrementFactor, bool anotation); //added incrementFactor 7-11-11

//R&D investment -- endogenous product innovation -- firms update current status by social learning with individual sampling of other firms
void updateRandDInvestmentStatus(firm firms[], const int firmNum, double intensityOfChoice, double probRandDInvestMutation, bool intermediateGoods, double overHeadLReq, double overHeadRandDReq, double salary, double intPriceDiscount, std::ofstream &testoutput, int round, bool anotation);

//R&D investment -- DC -- endogenous product innovation -- firms update current status by social learning with quasi global discrete choice
void updateRandDInvestmentStatusDC(firm firms[], const int firmNum, double intensityOfChoice, double weightRandDDC, double probRandDInvestMutation, bool intermediateGoods, double overHeadLReq, double overHeadRandDReq, double salary, double intPriceDiscount, double &medRecentRandDInvestment, double &avgRecentProfitAboveMedRandD, double &avgRecentProfitBelowMedRandD, double &avgRecentProfitType1, double &avgRecentProfitType2, std::ofstream &testoutput, int round, bool anotation);

//OH and R&D accounting -- //increment wageBill and firm profits for fixed overhead labor and current RandD overhead labor and increment recentRandDInvestment //10-7-10
void doOHandRandDAccounting(firm firms[], const int firmNum, double riGain, double overHeadLReq, double overHeadRandDReq, double salary, double &wageBill, double &salaryBill, int &numRandD, int &numRandD1, int &numRandD2, std::ofstream &testoutput, int round, bool anotation);

//labor productivity
double updateTechEff(firm firms[], const int firmNum, std::ofstream &testoutput, int &round, double &avTechEff, double &strtTech, bool anotation, double probEffShock);

//labor productivity -- common shock case
double updateTechEffCommonShock(firm firms[], const int firmNum, std::ofstream &testoutput, int &round, double &avTechEff, double &strtTech, bool anotation);

//markup pricing (from CES demands)
void updatePrice(firm firms[], const int firmNum, /*double inputFactor,*/ std::ofstream &testoutput, bool anotation, double intPriceDiscount);

//production and sale of intermediate goods
void doIntermediateGoods(firm firms[], double &wageBill, const int firmNum, double overHeadGReq, std::ofstream &testoutput, int round, bool anotation, double intPriceDiscount, bool intermediateGoods); 

//production of final (consumption) goods
void doProduction(firm firms[], double &wageBill, const int firmNum, double overHeadGReq, std::ofstream &testoutput, int round, bool anotation, double intPriceDiscount, bool intermediateGoods);

//utility -- calculate utility of representative consumer given a set of nominal firm shares in consumer spending
//changed wageBill to income to allow separation of wage and salary bill 5/13/16
double calcUtility(firm firms[], double income, const int firmNum, double shares[], std::ofstream &testoutput, int round, bool anotation, bool cesHedonics, double cesHedonicsElastGoods, double cesHedonicsElastChars, double complementIntConstraintFactor);

//consumer search -- set final goods demand shares
//added salaryBill 5/12/16
double doConsumerSearch(firm firms[], double wageBill, double salaryBill, const int firmNum, std::ofstream &testoutput, int round, int numMutationPairs, int numSearchTrials, double mutScale, bool anotation, bool cesHedonics, double cesHedonicsElastGoods, double cesHedonicsElastChars, double complementIntConstraintFactor);

//consumer search for two consumer/firm classes //5/13/16
double doConsumerSearchClass1(firm firms[], double wageBill, const int firmNum, std::ofstream &testoutput, int round, int numMutationPairs, int numSearchTrials, double mutScale, bool anotation, bool cesHedonics, double cesHedonicsElastGoods, double cesHedonicsElastChars, double complementIntConstraintFactor);

double doConsumerSearchClass2(firm firms[], double salaryBill, const int firmNum, std::ofstream &testoutput, int round, int numMutationPairs, int numSearchTrials, double mutScale, bool anotation, bool cesHedonics, double cesHedonicsElastGoods, double cesHedonicsElastChars, double complementIntConstraintFactor);

//final sales 	       
double doFinalSales(firm firms[], bool twoClasses, double wageBill, double salaryBill, const int firmNum, int firmNum1, int firmNum2, double rpGain, std::ofstream &testoutput, int round, bool anotation, bool altPriceIndex, double &nomGDP, double &realGDP, double &gdpDeflator, double &avgProfit, double &avgProfitType1, double &avgProfitType2); 

//restarts -- for hedonics 
int doRestartsHedonics(firm firms[], const int firmNum, double strtTech, double wageBill, double salaryBill, double overHeadGReq, double intPriceDiscount, bool imitation, double probRandomRestart, double probIntegerConstrained, double probComplements, double wage, double salary, double markup);

//restarts -- for hedonics with classes and with n1 n2 endogenous
int doRestartsHedonicsEndogN1N2(firm firms[], const int firmNum, int &firmNum1, int &firmNum2, double strtTech, double wageBill, double salaryBill, double overHeadGReq, double intPriceDiscount, bool imitation, double probRandomRestart, double probIntegerConstrained, double probComplements, double wage, double salary, double markup, double avgRecentProfitType1, double avgRecentProfitType2, double probSwitchMarkets, double probRandomSwitchMarkets, int ownMarketBias);


//restarts -- common shock case -- needs to be redone if used ***
int doRestartsCommonShock(firm firms[], const int firmNum, double strtTech, double wageBill, double overHeadGReq, double intPriceDiscount, bool imitation, double probRandomRestart, double wage, double salary, double markup); 

#endif
