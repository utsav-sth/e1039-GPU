# References

* https://github.com/E1039-Collaboration/e1039-core/blob/6f13d8621dd6577aa7427ba835ca12f28d04f5f6/framework/fun4all/Fun4AllDstInputManager.cc
 * https://github.com/E1039-Collaboration/e1039-core/blob/master/packages/reco/ktracker/Fun4AllSRawEventInputManager.cxx#L127
 * https://github.com/E1039-Collaboration/e1039-core/blob/6f13d8621dd6577aa7427ba835ca12f28d04f5f6/interface_main/SQSpillMap.h



# Data Structures

```cpp
class Spill {
	int SpillID;
	int RunID;
};

/**
 * @see https://github.com/E1039-Collaboration/e1039-core/blob/6f13d8621dd6577aa7427ba835ca12f28d04f5f6/interface_main/SQEvent.h
 */
class gSQEvent {
	public:
	int EventID;
	/**
	 * One spill can have many events.
	 */
	int SpillID;

	int RunId;

	int TriggerBits;
	short TargetPos;
	int TurnID;
	int RFID;
	int Intensity[33];
	short TriggerEmu;
	short NRoads[4];
	int NHits[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
	int nAH;
	int nTH;
	gHit AllHits[EstnAHMax];
	gHit TriggerHits[EstnTHMax];   
};

class SQEvent {
    // ...
    run_id() 
    spill_id() 
    event_id() 
    data_quality() 
    vme_time() 


    trigger(const SQEvent::TriggerMask i) 
    trigger() 
    raw_matrix(const unsigned short i) 
    after_inh_matrix(const unsigned short i) 
    qie_presum(const unsigned short i) 
    qie_trigger_count() 
    qie_turn_id() 
    qie_rf_id() 
    qie_rf_intensity(const short i) 
    flag_v1495() 
    n_board_qie() 
    n_board_v1495() 
    n_board_taiwan() 
    n_board_trig_bit() 
    n_board_trig_count() 

    // new properties
    coda_event_id()   // is unique to run (can get from runId)
}

// ######################################


class gHit {
	public:
	int index; // -> hit_id
	short detectorID; // 
	short elementID;
	float tdcTime;
	float driftDistance;
	float pos;

    // bool isInTime() const { return (flag & Hit::inTime) != 0; }
	short flag;
};



```

# Questions

* Get new ROOT file format.
* For the current online_reconstruction project, we don't seem to be using `LoadInputDict`.
    * Try: remove `LoadInputDict` from the makefile and see if it still works.
* Where are the following `Hit` properties?
    * `driftDistance`
    * `pos`
    * `flag`

## New Questions
* Is there a code sample that already reads that DST.root (that we have)?
    * Which of the `InputManagers` can read the file?
    * What is the file's `topnodename` etc?
    * How to add the libraries? (are there any pitfalls? where should I download them from etc?)

# Notes

* 
* All SQ* data types are implemented as PHObject
  * https://github.com/sPHENIX-Collaboration//coresoftware/blob/master/offline/framework/phool/PHObject.h


/**
 *
 */

https://github.com/sPHENIX-Collaboration//coresoftware/blob/master/offline/framework/phool/PHObject.h