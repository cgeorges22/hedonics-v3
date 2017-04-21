//Firm Class Declarations and Definitions

#ifndef FIRM_H
#define FIRM_H
#include <cassert> //7/16/09

const int numHedonicElements = 100; //10 //50 //100
const int maxNumComplements = 3; //added 6/26/09

class firm {
public:
  // Constructors
  firm() {
    active = 1;
    type = 1; //for segmented consumer market -- type 1 is goods for prod workers, type 2 is goods for OH labor // 5/11/16
    capital = 100.0; //default is 80 // has been 200 // for no restarts set very large
    //Qual = 1.0; // there is no longer an overall quality measure // removed 9/9/10 // may revive later
    techEff = 1.0; //1.0 // efficiency of production workers (marginal output per worker)
    effStdDev = 0.01; //0.01 // std dev of shocks to techEff
    techEffRandD = 1.0; // 1.0 // efficiency of RandD workers // not currently used -- ** should replace with overHeadlReq and overHeadRandDReq 
    inputFactor = 0.0; // 0.1 -- additional input needed from neighbor firm in units of neighbor's good per unit output // moved to firm.h 9/9/10
    wage = 1.0; //wage of production workers // currently not included in input.txt **
    //wageRandD = 1.0; //wage of RandD workers
    salary = 1.0; //salary of OH workers (general + RandD) // 5/12/16 // currently called from global not from here **
    price = 0.0; //price of firm's consumption good
    price0 = 0.0; //added 6/8/09 // 
    markup = 2.0; //moved here 9/9/10 // markup of price over MC //need markup>1 and markup*inputFactor*intPriceDiscount < 1
    //capDeficit = 0.0;
    //brrwCap = 0.0;
    finalDemand = 0.0;
    finalDemandShare = 0.0;
    exptdFinalDemand = 0.0;
    exptdIntermedDemand = 0.0;
    IntGoodOrder = 0.0;
    intGoodProduction = 0.0;
    finalProduction = 0.0;
    //invstmntVC = 0.0;
    //addProduction = 0.0;
    maxFinalProduction = 0.0;
    hasComplements = false; 
	numComplements = 0; 
    integerConstrained = false; 
	RandDInvestment = false;
	recentRandDInvestment = 0.0; 
	profit = 0.0;
	recentProfits = 0.0; 
    //hedonic arrays are initialized in main
 };

  //Get the (consumer good) type of the firm // 5/11/16
  int getType() const {return type;}

  //Get the capital of the firm
  double getCapital() const {return capital;}
  
  //Get the maximum the firm can produce
  double getMaxFinalProduction() const {return maxFinalProduction;}

  //Get the quality level of what the firm produces
  //double getQual() const {return Qual;}

  //Get the efficiency of production labor
  double getTechEff() const {return techEff;}

  //Get the standard deviation of the techEff inovation
  double getEffStdDev() const {return effStdDev;}
  
  //Get the fficiency of RandD labor
  double getTechEffRandD() const {return techEffRandD;}
  
  //Get the input factor for the firm (input from neighbor per unit output)
  double getInputFactor() const {return inputFactor;}

  //Get the wage of the firm's production labor
  double getWage() const {return wage;}
  
  //Get the wage of the firm's RandD labor
  //double getWageRandD() const {return wageRandD;}
  
  //Get the salary of the firm's OH labor
  double getSalary() const {return salary;}

  //Get the price the firm is charging
  double getPrice() const {return price;}
  
  //Get the initial price the firm charged in period 0 //added 6/8/09
  double getPrice0() const {return price0;}
  
  //Get the markup of price over MC //added 9/9/10
  double getMarkup() const {return markup;}

  //Get the amount of capital deficiency
  //double getCapDeficit() const {return capDeficit;}

  //Get the amount of capital borrowed // not used
  //double getBrrwCap() const {return brrwCap;}
  
  //Get the share of nominal demand for final goods //7/30/08
  double getFinalDemandShare() const {return finalDemandShare;}
  
  //Get the actual demand for the firm's consumption good
  double getFinalDemand() const {return finalDemand;}

  //Get the expected demand for the firm's consumption good //7/30/08
  double getExptdFinalDemand() const {return exptdFinalDemand;}

  //Get the expected demand firm's intermediate good
  double getExptdIntermedDemand() const {return exptdIntermedDemand;}

  //Get the order for the firm's intermediate good
  double getIntGoodOrder() const {return IntGoodOrder;}

  //Get the units of the intermediate good supplied
  double getIntGoodProduction() const {return intGoodProduction;}

  //Get the production of the consumption good
  double getFinalProduction() const {return finalProduction;}

  //Get the investment towards venture capital // not used
  //double getInvstmntVC() const {return invstmntVC;}

  //Get the market funded production // not used
  //double getAddProduction() const {return addProduction;}
  
  //Get the active status
  int getActive() const {return active;}
  
  //Get the hedonic quality of the firm's good //7/29/08
  int getHedonicQualityElementJ(int j) const {
	assert(j>=0 && j<numHedonicElements); //3/24/10
	return hedonicQuality[j];
  } 
  
  //Get is the firm's good indivisible in the consumption market // added 6/25/09
  bool getIntegerConstrained() const {return integerConstrained;} // is bool right here??
 
  //Get whether the firm has complements //6/26/09 //check bool here
  bool getHasComplements() const {return hasComplements;}

  //Get the number of the firm's complementary products //6/26/09
  int getNumComplements() const {return numComplements;}

  //Get complementary product //6/26/09
  int getComplementaryProductK(int k) const {
	assert(k>=0 && k< maxNumComplements); //7/16/09
	return complementaryProducts[k];
  }	

  //Get the hedonic quality from complementarity //6/26/09
  //int getComplementaryHedonicElementKJ(int k, int j) const {return complementaryHedonicElements[k][j];} //7/3/09  
  int getComplementaryHedonicElementKJ(int k, int j) const {
	assert(k>=0 && k<maxNumComplements); //7/16/09
	assert(j>=0 && j<numHedonicElements); //7/16/09
	assert(k * numHedonicElements + j >=0 && k * numHedonicElements + j < maxNumComplements * numHedonicElements); //7/16/09
	return complementaryHedonicElements[k * numHedonicElements + j]; //7/3/09
  }
  
  //Get whether the firm is currently innovating //3/20/10
  bool getRandDInvestment() const {return RandDInvestment;}
  
  //Get the firm's recent average incidence of innovation //3/20/10
  double getRecentRandDInvestment() const {return recentRandDInvestment;}
 
  //Get the firm's current profit //3/20/10
  double getProfit() const {return profit;}
      
  //Get the firm's recent average profits //3/20/10
  double getRecentProfits() const {return recentProfits;}

  //Set the market funded production // not used
  //void setAddProduction(double a) {addProduction = a;}
  
  //Set the firm's (consumer good) type // 5/11/16
  void setType(int a) {type = a;}

  //Set the maximum the firm can produce
  void setMaxFinalProduction(double a) {maxFinalProduction = a;}
  
  //Set the production of the consumption good
  void setFinalProduction(double a) {finalProduction = a;}

  //Set the investment in venture capital // not used
  //void setInvstmntVC(double a) {invstmntVC = a;}
  
  //Set the share of nominal demand for final goods //7/30/08
  void setFinalDemandShare(double a) {finalDemandShare = a;}
  
  //Set the actual final demand for consumption goods //7/30/08
  void setFinalDemand(double a) {finalDemand = a;}

  //Set the expected demand for the firm's consumption good
  void setExptdFinalDemand(double a) {exptdFinalDemand = a;}

  //Set the expected demand for the firm's intermediate good
  void setExptdIntermedDemand(double a) {exptdIntermedDemand = a;}

  //Set the order for the firm's intermediate good
  void setIntGoodOrder(double a) {IntGoodOrder = a;}

  //Set the firm's production of the intermediate good
  void setIntGoodProduction(double a) {intGoodProduction = a;}

  //Set the capital deficiency
  //void setCapDeficit(double a) {capDeficit = a;}

  //Set the capital borrowed // not used
  //void setBrrwCap(double a) {brrwCap = a;}

  //Set the capital of the firm
  void setCapital(double a) {capital = a;}

  //Set the quality level of what the firm produces
  //void setQual(double a) {Qual = a;}

  //Set the efficiency of production labor
  void setTechEff(double a) {techEff = a;}

  //Set the standard deviation of techEff inovation
  void setEffStdDev(double a) {effStdDev = a;}
  
  //Set the efficiency of RandD labor
  void setTechEffRandD(double a) {techEffRandD = a;}
  
  //Set the firm's input factor
  void setInputFactor(double a) {inputFactor = a;}

  //Set the wage of the firm's production labor
  void setWage(double a) {wage = a;}
  
  //Set the wage of the firm's production labor
  //void setWageRandD(double a) {wageRandD = a;}
  
  //Set the wage of the firm's production labor
  void setSalary(double a) {salary = a;}

  //Set the price the firm is charging
  void setPrice(double a) {price = a;}
  
  //Set the price the firm charges in initial round 0 //added 6/8/09
  void setPrice0(double a) {price0 = a;}
  
  //Set the markup of price over MC //9/9/10
  void setMarkup(double a) {markup = a;}

  //Set active status
  void setActive(int a) {active = 1;}
  
  //Set the hedonic quality of the firm's good //7/15/08
  void setHedonicQualityElementJ(int j, int a) {
	assert(j>=0 && j<numHedonicElements && a>=0); //3/24/10
	hedonicQuality[j] = a;
  } //7/15/08
  
  //Set is the firm integer constrained in the consumption market //added 6/25/09
  void setIntegerConstrained(bool a) {integerConstrained = a;} //is bool right here??

  //Set whether the firm has complements //6/26/09 //check bool here
  void setHasComplements(bool a) {hasComplements = a;}

  //Set the number of the firm's complementary products //6/26/98
  void setNumComplements(int a) {numComplements = a;}

  //Set complementary product //6/26/09
  void setComplementaryProductK(int k, int a) {
	assert(k>=0 && k< maxNumComplements && a>=0); //7/16/09
	complementaryProducts[k] = a;
  }	

  //Set the hedonic quality from complementarity //6/26/09
  //void setComplementaryHedonicElementKJ(int k, int j, int a) {complementaryHedonicElements[k][j] = a;}
  void setComplementaryHedonicElementKJ(int k, int j, int a) {
	assert(k>=0 && k<maxNumComplements); //7/16/09
	assert(j>=0 && j<numHedonicElements); //7/16/09
	assert(a>=0); //7/16/09
	assert(k * numHedonicElements + j >=0 && k * numHedonicElements + j < maxNumComplements * numHedonicElements); //7/16/09
	complementaryHedonicElements[k * numHedonicElements + j] = a;}
	
  //Set whether the firm is currently innovating //3/20/10
  void setRandDInvestment(bool a) {RandDInvestment = a;}

  //Set the firm's recent average incidence of innovation //3/20/10
  void setRecentRandDInvestment(double a) {recentRandDInvestment = a;}
  
  //Set the firm's current profit //3/20/10
  void setProfit(double a) {profit = a;}  
  
  //Set the firm's recent average profits //3/20/10
  void setRecentProfits(double a) {recentProfits = a;}

private:

  int active;
  int type; //5/11/16
  double maxFinalProduction;
  double IntGoodOrder;
  double capital;
  //double Qual;
  double techEff;
  double effStdDev;
  double techEffRandD;
  double inputFactor;
  double wage;
  //double wageRandD;
  double salary; //5/12/16
  double price;
  double price0;
  double markup;
  //double capDeficit;
  //double brrwCap;
  double finalDemandShare;
  double finalDemand;
  double exptdFinalDemand;
  double exptdIntermedDemand;
  double intGoodProduction;
  double finalProduction;
  //double invstmntVC;
  //double addProduction;
  int hedonicQuality[numHedonicElements]; //7/15/08 
  bool hasComplements; //6/26/09
  int numComplements;  //6/26/09
  int complementaryProducts[maxNumComplements]; //6/26/09
  //int complementaryHedonicElements[maxNumComplements][numHedonicElements]; //6/26/09
  int complementaryHedonicElements[maxNumComplements * numHedonicElements]; //6/26/09
  bool integerConstrained; //added 6/25/09
  bool RandDInvestment; //3/20/10
  double recentRandDInvestment; //3/20/10
  double recentProfits; //3/20/10
  double profit; //3/22/10
};

#endif
