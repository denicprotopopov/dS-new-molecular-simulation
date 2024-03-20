#pragma once

#include "ofMain.h"
#include "ofxGui.h"


#define INCLUDE_TRAILS false;

struct Particle {
    glm::vec2 position;
    glm::vec2 velocity;
    float mass;
    float kineticEnergy;
    float radius;
    std::deque<glm::vec2> trail;
    int trailMaxLength = 1;
    bool hasCollided;
};


class ofApp : public ofBaseApp{

	public:
		void setup() override;
		void update() override;
		void draw() override;
		void exit() override;
    
        ofxPanel gui;
        ofxFloatSlider targetTemperature;
        ofxFloatSlider coupling;
        ofxToggle applyThermostat;
        ofxIntSlider particleAmount;
        ofxButton applyAmount;
        ofxLabel label;
    

        vector<Particle> particles;
        float timeStep = 0.01;

        glm::vec2 computeForce(Particle &particle);
        void calculateTotalKineticEnergy(Particle &particle);
        void updateParticle(Particle &particle, float deltaTime);
        void checkWallCollisions(Particle &particle);
        void applyBerendsenThermostat();
        void resolveParticleCollisions();
        void initializeParticles();
		void keyPressed(int key) override;
};
