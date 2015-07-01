#!/usr/bin/ruby
include Process

$exe = "a.out"
$maxRuns = 8
$running = 0
$rCount = 0

thickness = [100, 200, 300, 400, 500, 600]
convective = [0, 1]
densityGradient = [3e-7, 4e-7, 5e-7]

runs = thickness.product(convective, densityGradient)

def start(d, g, c)
	cmd = "./" + $exe + " #{d} #{g} #{c}"
	pid = fork do
		exec cmd
	end
	puts "(run=%03d) %s started (pid=%05d)" % [$rCount, cmd, pid]
	return pid
end

pids = []

while 1 do
	if (pids.length < $maxRuns)
		x = runs[$rCount][0]
		y = runs[$rCount][1]
		z = runs[$rCount][2]
		pids << start(x, z, y)
		$rCount += 1
	else
		pid = wait()
		puts ""
		puts "%s -- %05d finished" % [Time.now, pid]
		pids.delete_if{|p| p == pid}
	end
end
