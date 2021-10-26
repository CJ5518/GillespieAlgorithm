function restoreState()
	state = {20, 0, 480}
	infectRate = .5;
	recoveryRate = .1;
	time = 0;
end

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

stoichiometry = {
	{1, 0, -1}, --sus to infected
	{-1, 1, 0} -- infected to recoved
}
--Infected, recoved, susd
--https://github.com/wefatherley/monte-carlo
state = {}
function normalTick(dt)
	print("dt", dt)
	local newState = {};
	for i, v in pairs(state) do
		newState[i] = v;
	end
	for i, v in pairs(propensity) do
		local number = v();
		for i2, v2 in pairs(stoichiometry[i]) do
			--print(i, i2, number, v2);
			newState[i2] = newState[i2] + (number * dt * v2);
		end
	end
	state = newState;
	time = time + dt;

	return true;
end

--https://en.wikipedia.org/wiki/Poisson_distribution#Generating_Poisson-distributed_random_variables
function poissonNumber(lambda)
	local L = math.exp(-lambda);
	local k = 0;
	local p = 1;
end

--Doesn't do the propensity cylcing thingy
function gillespieTick()
	--Summing everything gives us a result of 0, so something is messed up majorly
	local totalSum = 0;
	for i, v in pairs(propensity) do
		totalSum = totalSum + v();
	end
	--Nothing will happen
	if totalSum == 0 then return false end;
	local otherThing = math.log(1 / math.random());
	local sojourn = otherThing / totalSum;
	
	local thingToBeat = totalSum * math.random();

	local index = 0;
	while thingToBeat >= 0 do
		index = index + 1;
		thingToBeat = thingToBeat - propensity[index]();
	end

	for i, v in pairs(stoichiometry[index]) do
		state[i] = state[i] + v;
	end
	time = time + sojourn;
	return true;
end

--Assume this does something
function possionNumber(a, tau)
	return 1;
end

--https://aip.scitation.org/doi/pdf/10.1063/1.1378322
function tauLeaping()

end

function logRun(runs, logFunc, tickFunc, ...)
	local sus,infect, recov, dts = {},{},{},{};
	restoreState();
	local function logState()
		sus[#sus+1] = state[3];
		infect[#infect+1] = state[1];
		recov[#recov+1] = state[2];
		dts[#dts+1] = time;
	end
	for q = 1, runs do
		logState()
		if not tickFunc(...) then
			break;
		end
	end
	logFunc(sus,infect,recov,dts);
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

function logToCSV(sus,infect, recov, dts)
	local file = io.open("Output.csv", "w");
	file:write("dt,Susceptible,Infected,Recovered\n");
	for i, v in ipairs(dts) do
		local str = string.format("%f,%f,%f,%f\n", dts[i], sus[i], infect[i], recov[i]);
		file:write(str);
	end
	file:close();
end

restoreState();
logRun(1000, logToCSV, gillespieTick)