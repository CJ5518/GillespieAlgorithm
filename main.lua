require("gamma")

--[[
Examples:
logRun(7, logToConsole, {}, gillespieTick, {})
logRun(500, logToCSV, {}, tauLeaping, {.1, nil})
logRun(500, logToCSV, {}, tauLeaping, {nil, .1})
]]

--TODO:
--Better choose the value for tau

math.randomseed(os.time());
--Tis an old wives tale
math.random();
math.random();
math.random();

propensity = {
	--Susceptible to infected
	function()
		return (infectRate * state[3] * state[1]) / (state[1] + state[2] + state[3]);
	end,
	--infected to recovered
	function()
		return recoveryRate * state[1];
	end,
}
--stoichiometry[reaction][state]
stoichiometry = {
	{1, 0, -1}, --sus to infected
	{-1, 1, 0} -- infected to recoved
}
--Infected, recoved, susd
--https://github.com/wefatherley/monte-carlo
state = {}

function restoreState()
	state = {200, 0, 2e6 - 200}
	infectRate = .5;
	recoveryRate = .1;
	time = 0;
end


function cloneState()
	local newState = {};
	for i, v in pairs(state) do
		newState[i] = v;
	end
	return newState;
end

--Basic dumb method with floating point people
--args[1] is dt
function normalTick(args)
	local dt = args[1] or 1;
	local newState = cloneState();
	for i, v in pairs(propensity) do
		local number = v();
		for i2, v2 in pairs(stoichiometry[i]) do
			newState[i2] = newState[i2] + (number * dt * v2);
		end
	end
	state = newState;
	time = time + dt;

	return true;
end

--The probability that k is the result given lambda
function poissonProbability(k, lambda)
	return math.exp((k * math.log(lambda)) - lambda - log_gamma(k+1));
end

--Pick your poisson
--O(n) time, but I didn't feel like implementing a faster method
function poissonNumber(lambda)
	local L = math.exp(-lambda);
	local p = 1;
	local k = 0;
	repeat
		k = k + 1;
		p = p * math.random();
	until not (p > L);
	return k - 1;
end

--Doesn't do the propensity cylcing thingy
--Which may or may not be necessary
--Basic gillespie algorithm
function gillespieTick()
	--Summing everything gives us a result of 0, so something is messed up majorly
	local totalSum = 0;
	for i, v in pairs(propensity) do
		totalSum = totalSum + v();
	end
	--Nothing will happen
	if totalSum == 0 then return false end;

	--Whens the reaction?
	local otherThing = math.log(1 / math.random());
	local sojourn = otherThing / totalSum;
	
	--Which reaction?
	local thingToBeat = totalSum * math.random();
	local index = 0;
	while thingToBeat >= 0 do
		index = index + 1;
		thingToBeat = thingToBeat - propensity[index]();
	end
	--Update state
	for i, v in pairs(stoichiometry[index]) do
		state[i] = state[i] + v;
	end
	time = time + sojourn;
	return true;
end

--https://aip.scitation.org/doi/pdf/10.1063/1.2159468
function chooseTauFancy(epsilon)
	epsilon = epsilon or .1;
	--Returns the highest order reaction the given state is a part of
	--I think I'm calculating this right?
	local function HOR(stateIdx)
		--0 order reaction lmao
		local order = 0;

		for i, v in pairs(stoichiometry) do
			if v[stateIdx] ~= 0 then order = order + 1; end
		end
		
		return order;
	end
	--The two hat funcs, encapsulated into one
	--Set sigma to true to use the sigma function
	local function hatFunc(stateIdx, isSigma)
		local sum = 0;
		for i, v in pairs(propensity) do
			local stoich = stoichiometry[i][stateIdx];
			--Pretty much just absolute value in our case, but the paper says square it
			if isSigma then stoich = stoich * stoich end
			sum = sum + (stoich * v());
		end
		return sum;
	end

	--Pray that we understood the HOR
	local function gFunc(reactionIdx)
		local HORRes = HOR(reactionIdx);
		if HORRes == 1 then 
			return 1;
		elseif HORRes == 2 then
			return 2;
		end
	end

	local tau = 1000000;

	for i, v in pairs(state) do
		tau = math.min(tau,
			math.min(
				math.max((epsilon * v) / gFunc(i), 1) / math.abs(hatFunc(i, false)), math.pow(math.max((epsilon * v) / gFunc(i), 1),2) / hatFunc(i, true)
			)
		)
	end
	return tau;
end

--https://aip.scitation.org/doi/pdf/10.1063/1.1378322
--args[1]:number - should be a fixed time to jump by, or nil if we should pick it ourself
--args[2]:number - (0,1) epsilon, the error in our tau leaping, if we pick it
function tauLeaping(args)
	--Select the value for tau
	local tau = args[1] or chooseTauFancy(args[2]);

	local newState = cloneState();

	local allWereZero = true;

	for i, v in pairs(propensity) do
		local propensityResult = v();
		local reactionEvents = poissonNumber(propensityResult * tau);

		if propensityResult > 0 then allWereZero = false; end

		for i2, v2 in pairs(stoichiometry[i]) do
			newState[i2] = newState[i2] + (reactionEvents * v2);
		end
	end
	state = newState;
	time = time + tau;
	return not allWereZero;
end

--args[1]:number - should be a fixed time to jump by, or nil if we should pick it ourself
--args[2]:number - (0,1) epsilon, the error in our tau leaping, if we pick it
function estimatedMidpointTauLeaping(args)
	local tau = args[1] or chooseTauFancy(args[2]);
	local expectedStateChange = {0,0,0}

	--Calculate the expected state
	for i, v in pairs(propensity) do
		local aj = v();
		for i2, v2 in pairs(stoichiometry[i]) do
			expectedStateChange[i2] = expectedStateChange[i2] + (aj * v2);
		end
	end

	local midpointState = {0,0,0};
	for i, v in pairs(expectedStateChange) do
		midpointState[i] = state[i] + (v / 2);
	end

	newState = cloneState();
	state = midpointState;

	local allWereZero = true;

	--Calculate the actual state
	for i, v in pairs(propensity) do
		local propensityResult = v();
		local reactionEvents = poissonNumber(propensityResult * tau);

		if propensityResult > 0 then allWereZero = false; end

		for i2, v2 in pairs(stoichiometry[i]) do
			newState[i2] = newState[i2] + (reactionEvents * v2);
		end
	end
	state = newState;
	time = time + tau;
	return not allWereZero;
end

function logRun(runs, logFunc, logFuncArgs, tickFunc, tickFuncArgs)
	local sus,infect, recov, dts = {},{},{},{};
	restoreState();
	--Log the state into the tables
	local function logState()
		sus[#sus+1] = state[3];
		infect[#infect+1] = state[1];
		recov[#recov+1] = state[2];
		dts[#dts+1] = time;
	end
	--Briefly check if things have gone negative
	local function isStateGood()
		for i,v in pairs(state) do
			--Set the negative things to zero
			if v < 0 then state[i] = 0 end
		end
		return true;
	end
	for q = 1, runs do
		logState()
		if not tickFunc(tickFuncArgs) then
			break;
		end
		if not isStateGood() then break end
	end
	logState();
	logFunc(sus,infect,recov,dts, logFuncArgs);
end

function logToConsole(sus,infect, recov, dts)
	print("-------------SUS--------------");
	for i, v in pairs(sus) do
		print(v);
	end
	print("-------------INFECTED--------------");
	for i, v in pairs(infect) do
		print(v);
	end
	print("-------------RECOVERED--------------");
	for i, v in pairs(recov) do
		print(v);
	end
	print("-------------DT--------------");
	for i, v in pairs(dts) do
		print(v);
	end
end

--args[1] should be the filename of the output file
function logToCSV(sus,infect, recov, dts, args)
	local file = io.open(args[1] or "Output.csv", "w");
	file:write("dt,Susceptible,Infected,Recovered\n");
	for i, v in ipairs(dts) do
		local str = string.format("%f,%f,%f,%f\n", dts[i], sus[i], infect[i], recov[i]);
		file:write(str);
	end
	file:close();
end

restoreState();
--logRun(1000, logToCSV, gillespieTick)