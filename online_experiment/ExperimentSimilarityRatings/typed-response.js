jsPsych.plugins["typed-response"] = (function() {

    var plugin = {};
  
    plugin.info = {
      name: "typed-response",
      parameters: {
        prompt: {
          type: jsPsych.plugins.parameterType.STRING,
          default: ""
        },
        trial_duration: {
          type: jsPsych.plugins.parameterType.INT,
          default: null,
          description: "The maximum amount of time to wait for a response."
        }
      }
    }
  
    plugin.trial = function(display_element, trial) {
  
      // create html for prompt and response
      var html = "<p>" + trial.prompt + "</p>";
      html += "<input type='text' id='response' autofocus>";
  
      // display prompt and response
      display_element.innerHTML = html;
  
      // set up variables for tracking response and timer
      var response = "";
      var timer;
  
      // start timer if trial duration is specified
      if (trial.trial_duration !== null) {
        timer = setTimeout(end_trial, trial.trial_duration);
      }
  
      // function to end the trial and save the response
      function end_trial() {
        // clear timer
        clearTimeout(timer);
  
        // save response
        response = document.getElementById("response").value;
  
        // data to be saved with the trial
        var trial_data = {
          response: response
        };
  
        // end trial and send data to jsPsych
        jsPsych.finishTrial(trial_data);
      }
  
      // add event listener to response input to end trial on enter key press
      document.getElementById("response").addEventListener("keydown", function(e) {
        if (e.key === "Enter") {
          end_trial();
        }
      });
  
      // select response input as soon as the trial starts
      document.getElementById("response").focus();
  
    };
  
    return plugin;
  
  })();