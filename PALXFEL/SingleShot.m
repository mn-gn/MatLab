classdef SingleShot

    properties
        TimeStamp   (1,1) double
        PID         (1,1) int64
        AI_int      (:,1) double
        I0          (1,1) double
        pumpOn      (1,1) logical
        AUC         (1,1) double
        NormFactor  (1,1) double

    end

    methods
        function obj = SingleShot(TimeStamp,PID,AI_int,is_on)
            if nargin > 0
                obj.TimeStamp   = TimeStamp;
                obj.PID         = PID;
                obj.AI_int      = AI_int;
                obj.pumpOn      = is_on;
            end
        end

     
    end
end
