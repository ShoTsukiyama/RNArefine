# Setting of compile
CXX = g++
#CXX = FCCpx
#CXX = mpiFCCpx
#CXXFLAGS = -O3 -std=c++17
#CXXFLAGS = -Nclang -std=c++17 -g -pthread
CXXFLAGS = -std=c++17 -g -pthread

OUTPUT_DIR = ../bin
TARGET = $(OUTPUT_DIR)/rnaoptimizer
#SRC_FILES = Main.cpp SSPred.cpp ParseSeq.cpp RNArefine.cpp LoadPDB.cpp LoadBaseTemplate.cpp CommonPara.cpp LoadParams.cpp BasicFunc.cpp SimulationBase.cpp SimulationCalcParams.cpp SimulationCalcCoord.cpp SimulationSaveStr.cpp SimulationEnergyMain.cpp SimulationEnergyFAPE.cpp SimulationEnergyClash.cpp SimulationEnergyBondLength.cpp SimulationEnergyAdjBondLength.cpp SimulationEnergyBondAngle.cpp SimulationEnergyTorsionAngle.cpp SimulationEnergyBasePair.cpp CalcBaseTemplateRMSD.cpp SimulationMovementMain.cpp SimulationMovementTrosion.cpp SimulationMovementBondAngle.cpp SimulationMovementBondLength.cpp SimulationMovementBaseCoord.cpp SimulationMovementLMP.cpp SimulationMovementStemCoord.cpp
SRC_FILES = CommonPara.cpp Main.cpp ParseSeq.cpp RNArefine.cpp LoadPDB.cpp LoadParams.cpp LoadDeepParams.cpp BasicFunc.cpp SimulationBase.cpp SimulationCalcParams.cpp SimulationCalcCoord.cpp SimulationSaveStr.cpp SimulationEnergyMain.cpp SimulationEnergyTemplate.cpp SimulationEnergyContact.cpp SimulationEnergyClash.cpp SimulationEnergyBondLength.cpp SimulationEnergyBondAngle.cpp SimulationEnergyTorsionAngle.cpp SimulationEnergyPairTorsionAngle.cpp SimulationEnergyPairBBTorsionAngle.cpp SimulationEnergyMultiTorsionAngle.cpp SimulationEnergyTorsionAngleCls.cpp SimulationEnergyBasePair.cpp SimulationEnergyBasePairHyd.cpp SimulationMovementGradHydGeometries.cpp SimulationEnergyStacking.cpp SimulationEnergySingleNtTemplate.cpp SimulationEnergyUnonCBasePairTemplate.cpp SimulationEnergySingleStackTemplate.cpp BaseTrans.cpp SimulationMovementMain.cpp SimulationMovementTorsion.cpp SimulationMovementBondAngle.cpp SimulationMovementBondAngle2.cpp SimulationMovementBondAngleMt.cpp SimulationMovementBondLength.cpp SimulationMovementBondLength2.cpp SimulationMovementBondLengthMt.cpp SimulationMovementBondLengthAngle.cpp SimulationMovementTorsionSugarRing.cpp SimulationMovementMultiTorsion.cpp SimulationMovementMultiTorsionSugar.cpp SimulationMovementRotate.cpp SimulationMovementLMP.cpp SimulationMovementLMPfull.cpp SimulationMovementLMPfullbasepair.cpp SimulationMovementgradientLength2.cpp SimulationMovementgradientLength3.cpp SimulationMovementgradientLength4.cpp GSSfunc.cpp SimulationMovementgradientAngle.cpp SimulationMovementHydGeometries.cpp Gradfunc.cpp SimulationMovementGeometries.cpp StrVis.cpp

OBJ_DIR = ./obj
OBJ_FILES = $(patsubst %.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
OBJ_NODIR_FILES = $(SRC_FILES:.cpp=.o)

CFLAGS += $(shell pkg-config --cflags eigen3)
LDLIBS += $(shell pkg-config --libs eigen3)
#LDLIBS += -lstdc++fs
#LDLIBS += -lc++fs

#all:
#echo $(OBJ_FILES)

$(TARGET):  $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(CFLAGS) -o $@ $^  $(LDLIBS)

$(OBJ_DIR)/CommonPara.o: CommonPara.cpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/RNArefine.o: RNArefine.cpp RNArefine.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/LoadPDB.o: LoadPDB.cpp LoadPDB.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/LoadParams.o: LoadParams.cpp LoadParams.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/LoadDeepParams.o: LoadDeepParams.cpp LoadDeepParams.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/ParseSeq.o: ParseSeq.cpp ParseSeq.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/BasicFunc.o: BasicFunc.cpp BasicFunc.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationBase.o: SimulationBase.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationCalcParams.o: SimulationCalcParams.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationCalcCoord.o: SimulationCalcCoord.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationSaveStr.o: SimulationSaveStr.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyTemplate.o: SimulationEnergyTemplate.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyContact.o: SimulationEnergyContact.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyClash.o: SimulationEnergyClash.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyBondLength.o: SimulationEnergyBondLength.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyBondAngle.o: SimulationEnergyBondAngle.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyTorsionAngle.o: SimulationEnergyTorsionAngle.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyPairTorsionAngle.o: SimulationEnergyPairTorsionAngle.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyPairBBTorsionAngle.o: SimulationEnergyPairBBTorsionAngle.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyMultiTorsionAngle.o: SimulationEnergyMultiTorsionAngle.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyTorsionAngleCls.o: SimulationEnergyTorsionAngleCls.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

#$(OBJ_DIR)/SimulationEnergyNtpairlength.o: SimulationEnergyNtpairlength.cpp SimulationBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

#$(OBJ_DIR)/SimulationEnergyAdjpairlength.o: SimulationEnergyAdjpairlength.cpp SimulationBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyBasePair.o: SimulationEnergyBasePair.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyBasePairHyd.o: SimulationEnergyBasePairHyd.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementGradHydGeometries.o: SimulationMovementGradHydGeometries.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyStacking.o: SimulationEnergyStacking.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergySingleNtTemplate.o: SimulationEnergySingleNtTemplate.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyUnonCBasePairTemplate.o: SimulationEnergyUnonCBasePairTemplate.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergySingleStackTemplate.o: SimulationEnergySingleStackTemplate.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

#$(OBJ_DIR)/SimulationEnergyPucker.o: SimulationEnergyPucker.cpp SimulationBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/BaseTrans.o: BaseTrans.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationEnergyMain.o: SimulationEnergyMain.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementMain.o: SimulationMovementMain.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementTorsion.o: SimulationMovementTorsion.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementBondAngle.o: SimulationMovementBondAngle.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementBondAngle2.o: SimulationMovementBondAngle2.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementBondAngleMt.o: SimulationMovementBondAngleMt.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementBondLength.o: SimulationMovementBondLength.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementBondLength2.o: SimulationMovementBondLength2.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementBondLengthMt.o: SimulationMovementBondLengthMt.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementBondLengthAngle.o: SimulationMovementBondLengthAngle.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementTorsionSugarRing.o: SimulationMovementTorsionSugarRing.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementMultiTorsion.o: SimulationMovementMultiTorsion.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementMultiTorsionSugar.o: SimulationMovementMultiTorsionSugar.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementRotate.o: SimulationMovementRotate.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementLMP.o: SimulationMovementLMP.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementLMPfull.o: SimulationMovementLMPfull.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementLMPfullbasepair.o: SimulationMovementLMPfullbasepair.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

#$(OBJ_DIR)/SimulationMovementgradientLength.o: SimulationMovementgradientLength.cpp SimulationBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementgradientLength2.o: SimulationMovementgradientLength2.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementgradientLength3.o: SimulationMovementgradientLength3.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementgradientLength4.o: SimulationMovementgradientLength4.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/GSSfunc.o: GSSfunc.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementgradientAngle.o: SimulationMovementgradientAngle.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/SimulationMovementHydGeometries.o: SimulationMovementHydGeometries.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/Gradfunc.o: Gradfunc.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@
	
$(OBJ_DIR)/SimulationMovementGeometries.o: SimulationMovementGeometries.cpp SimulationBase.hpp CommonPara.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/StrVis.o: StrVis.cpp StrVis.hpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

$(OBJ_DIR)/Main.o: Main.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJ_FILES)

#$(OBJ_DIR)/SimulationMovementStemCoord.o: SimulationMovementStemCoord.cpp SimulationBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

#$(OBJ_DIR)/CalcBaseTemplateRMSD.o: CalcBaseTemplateRMSD.cpp SimulationBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) $(CFLAGS) $(LDLIBS) -c $< -o $@

#$(OBJ_DIR)/SimulationMovementMain.o: SimulationMovementMain.cpp SimulationBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalBase.o: SimulationLocalBase.cpp SimulationLocalBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalEnergyBondAngle.o: SimulationLocalEnergyBondAngle.cpp SimulationLocalBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalEnergyTorsionAngle.o: SimulationLocalEnergyTorsionAngle.cpp SimulationLocalBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalEnergyBondLength.o: SimulationLocalEnergyBondLength.cpp SimulationLocalBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalEnergyClashHbond.o: SimulationLocalEnergyClashHbond.cpp SimulationLocalBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalEnergyMain.o: SimulationLocalEnergyMain.cpp SimulationLocalBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalCalcParams.o: SimulationLocalCalcParams.cpp SimulationLocalBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalCalcCoord.o: SimulationLocalCalcCoord.cpp SimulationLocalBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalMovementMain.o: SimulationLocalMovementMain.cpp SimulationLocalBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalMovementTrosion.o: SimulationLocalMovementTrosion.cpp SimulationLocalBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalMovementBondAngle.o: SimulationLocalMovementBondAngle.cpp SimulationLocalBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalMovementBondLength.o: SimulationLocalMovementBondLength.cpp SimulationLocalBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalSaveStr.o: SimulationLocalSaveStr.cpp SimulationLocalBase.hpp CommonPara.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJ_DIR)/SimulationLocalEnergyClashHbond.o: SimulationLocalEnergyClashHbond.cpp SimulationLocalBase.hpp
#@mkdir -p $(OBJ_DIR)
#$(CXX) $(CXXFLAGS) -c $< -o $@

#$(TARGET): $(OBJ_FILES)
#$(CXX) $(CXXFLAGS) $^ -o $@

# Build rules for each CPP file to generate object files
#%.o: %.cpp
#@mkdir -p $(@D)
#$(CXX) $(CXXFLAGS) -c $< -o $(OBJ_DIR)/$@

#link: $(OBJ_FILES)
#@$(CXX) $(CXXFLAGS) $< -o $(TARGET)

# Custom link target to link the object files

#.PHONY: all link clean

