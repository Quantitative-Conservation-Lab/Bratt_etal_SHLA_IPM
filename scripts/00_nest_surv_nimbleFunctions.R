# AEB
# Sept 15, 2020
# exploring nimble functions to write stage-based nest survival model

#### GET STATE ####
GetStateTimeVary <- nimbleFunction(
  run = function(state = double(0),
                 phi = double(0),
                 laymin =  double(0),
                 laymax = double(0),
                 incmin = double(0),
                 incmax = double(0),
                 nesmin = double(0),
                 nesmax = double(0),
                 timeInCurrState = double(0),
                 etaLay = double(1),
                 etaInc = double(1),
                 etaNest = double(1)
  ) 
  {
    
    returnType(double(1)) # returns appropriate row vector
    ps <- rep(0, 5)   # set default as zero, then just fill in non-zero cells
    
    # STATES (alive)
    # 1 - laying
    # 2 - incubating
    # 3 - nestling
    # 4 - fledged
    # 5 - failed
    
    if (state == 1 & timeInCurrState < laymax) { # if currently laying and could still be laying
      ps[1] <- phi * (etaLay[timeInCurrState])
      ps[2] <- phi * (1 - etaLay[timeInCurrState])
      ps[5] <- (1 - phi)
    } 
    
    if (state == 1 & timeInCurrState >= laymax) { # if already laying and must transition
      ps[1] <- 0
      ps[2] <- phi 
      ps[5] <- (1 - phi)
    } 
    
    if (state == 2 & timeInCurrState < incmin) { # if currently incubating and MUST still be
      ps[2] <- phi
      ps[3] <- 0
      ps[5] <- (1 - phi)
    }
    
    if (state == 2 & timeInCurrState >= incmin & timeInCurrState < incmax) { # if currently incubating and COULD still be
      ps[2] <- phi * (etaInc[timeInCurrState - incmin + 1]) # AEB - to do - change this in the model function
      ps[3] <- phi * (1 - etaInc[timeInCurrState - incmin + 1])
      ps[5] <- (1 - phi)
    }
    
    if (state == 2 & timeInCurrState >= incmax) { # if currently incubating and must transition
      ps[2] <- 0
      ps[3] <- phi
      ps[5] <- (1 - phi)
    }
    
    if (state == 3 & timeInCurrState < nesmin) { # if currently nestling and MUST still be
      ps[3] <- phi
      ps[4] <- 0
      ps[5] <- (1 - phi)
    }
    
    if (state == 3 & timeInCurrState >= nesmin & timeInCurrState < nesmax) { # if currently nestling and COULD still be
      ps[3] <- phi * (etaNest[timeInCurrState - nesmin + 1])
      ps[4] <- phi * (1 - etaNest[timeInCurrState - nesmin + 1])
      ps[5] <- (1 - phi)
    }
    
    if (state == 3 & timeInCurrState >= nesmax) { # if currently nestling and must transition
      ps[3] <- 0
      ps[4] <- phi
      ps[5] <- (1 - phi)
    }
    
    if (state == 4) { # fledged stay fledged
      ps[4] <- 1
    }
    
    if (state == 5) { # dead stay dead
      ps[5] <- 1
    }
    
    return(ps[1:5])
  }
)

CGetStateTimeVary <- compileNimble(GetStateTimeVary)

#### GET TIME SINCE LAST TRANSITION ####
GetTimeSinceLastTransition <- nimbleFunction(
  run = function(state=double(1), 
                 first=double(0),
                 firsttime = double(0),
                 d=double(0)
  ) 
  {
    returnType(double(0))
    
    ans <- 0
    if (state[d-first+1] == state[1]) { # if has not yet transitioned
      ans <- firsttime + d - first
    } else { # if state has transitioned
      # TODO do we need to add one here?
      ans <- length(which(state[1:(d-first+1)] == state[d-first+1]))
    }
    
    return(ans)
  }
)

#### GET TIME IN STATE AT INITIAL DETECTION ####
GetTimeinStateatInitDetection <- nimbleFunction(
  run = function(state=double(0), 
                 TimeInStageInitDetectionLay=double(0),
                 TimeInStageInitDetectionInc=double(0),
                 TimeInStageInitDetectionNes=double(0)
  ) 
  {
    returnType(double(0))
    
    ans <- 0
    if (state == 1) {
      ans <- TimeInStageInitDetectionLay
    } else if (state == 2) {
      ans <- TimeInStageInitDetectionInc
    } else if (state == 3) {
      ans <- TimeInStageInitDetectionNes
    }
    
    return(ans)
  }
)

CGetTimeSinceLastTransition <- compileNimble(GetTimeSinceLastTransition)