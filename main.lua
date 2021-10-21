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
	local newState = {};
	for i, v in pairs(propensity) do
		newState[i] = state[i] + (v() * dt);
	end
	state = newState;
	time = time + dt;
end

function gillespieTick()
	--Summing everything gives us a result of 0, so something is messed up majorly
	local totalSum = propensity[2]() - propensity[3]();
	local otherThing = math.log(1 / math.random());
	local sojourn = otherThing / totalSum;
	
	--The j part



end

function printState()
	print(susceptible, infected, recovered);
end

function logRun(step, runs)
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
		normalTick(step);
	end

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

restoreState();