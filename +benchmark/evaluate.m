function result = evaluate(prediction)
%EVALUATE Third module: participant-level clinical/control benchmark contract.
result=tbench.benchmark.metrics(prediction);
end
