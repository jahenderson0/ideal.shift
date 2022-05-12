# ideal.shift
############################################
Some prototype code for measuring issue-by-issue ideal points of MCs in the 110th Congress

Name: John Henderson
Data: April 18, 2019 [updated May 11, 2022]
############################################


GOAL: produce a model of ideal points that can shift isssue-by-issue; reflects variation in how MCs represent districts depending on the issue areas at stake
APPROACH: scale cosponsor choices with top-coded issue areas using the Policy Agendas Project codings; see here: http://congressionalbills.org/index.html
  - to replicate topic mixtures, the idea is use words in bill titles to predict PAP issue areas, this gives topic proportion estimates for each bill
  - these topic proportions are used as weights in an unsupervised IRT model, so that ideal points are allowed to float left or right on cosponsorship choices that have more topic weight given to topic
  - the cutpoint model then is Pr(Cosponsor=1) = f{(alpha+t(z)%*%phi)*beta + gamma + delta};
  - where alpha are ideal points over all choices; z reflect topic-specific shifts in ideal points 
  - beta are bill discrimination parameters; delta are bill location parameters; gamma are legislator connectedness paramaters
  - a threshold model is also explored: Pr(Cosponsor=1) = f{(alpha+t(z)%*%phi - beta)^2 + gamma + delta}
